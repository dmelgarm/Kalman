function [dk vk]=kalmans(d,a,Td,Ta,q,r)

% D.Melgar 07/2010
% 
% Perform Kalman filter on accelerometer and GPS data to obtain 'true' 
% displacements and velocities (see Smyth & Wu, 2007). The time series have
% to be aligned on the FIRST sample and da/dd must be an integer. The
% acceleration time series can be longer than the gps time series because 
% position will jsut be computed with time update but NOT viceversa.
% This routine includes a Rauch-Trung-Striebel smoothing algorithm
% 
% IN
% d - displacement time series
% a - acceleration time series
% Td -  displacement sample rate in seconds
% Ta - acceleration sample rate in seconds
% q - process noise variance (Accel.)
% r -  measurement noise variance (GPS)
%
% 
% OUT
% dk - kalman filtered displacement
% vk  - kalman filtered velocity

% DMM 10/2010
%
% Added memory preallocation to speedup runtime you need functions rnd() 
% and isintdiv() to execute.


%Structural preeliminaries
if size(d,1)<size(d,2)  %turn to column vector
    d=d';
end
if size(a,1)<size(a,2)  %ditto
    a=a';
end
la=size(a,1);
sizegps=size(d,1);
%Initalize output vectors
dk=zeros(la,1);
vk=zeros(la,1);
xplus=zeros(2,la);
xminus=zeros(2,la);
s=zeros(2,la);
Pplus=zeros(2,2,la);
Pminus=zeros(2,2,la);


%Initialize state variable matrices 
x=[0;0];
%Initialize observation matrix (displacement)
z=0;
%Initialize state covariance matrix
P=[0 0;0 0];
%Initialize Kalman gain matrix
K=[0 0;0 0];
%Initialize State-Space matrices
A=[1 Ta;0 1];
B=[Ta^2/2;Ta];
H=[1 0];
%Initalize noise covariance matrices
%Q=[q*Ta^3/3 q*Ta^2/2;q*Ta^2/2 q*Ta];
%R=r/Td;
%Acceleration to dispalcement sampling ratio
ratio=rnd(Td/Ta);
if max(size(q))==1
    q=repmat(q,1,la);
end
if max(size(r))==1
    r=repmat(r,1,sizegps);
end

%Kalmanize!

%FORWARD FILTER
%k counts acceleration data, i counts gps data
i=1;
measure=0;
for k=1:la
    Q=[q(k)*Ta^3/3 q(k)*Ta^2/2;q(k)*Ta^2/2 q(k)*Ta];    %TRIAL CODE
    %Time update
    %predict state
    x=A*x+B*a(k);
    %predict covariance
    P=A*P*A'+Q;
    xminus(:,k)=x;
    Pminus(:,:,k)=P;   
    %Measurement update (if GPS is available)
    if i< sizegps & isintdiv(k-1,ratio)==1 & ~isnan(d(i))
       R=r(i)/Td;    %TRIAL CODE
       %compute Kalman gain
       K=P*H'/(H*P*H'+R);
       %read in GPS displacement
       %update state
       x=x+K*(d(i)-H*x);
       %update covariance
       %P=(eye(2)-K*H)*P;
       P=(eye(2)-K*H)*P*(eye(2)-K*H)'+K*R*K';
       i=i+1;
       Pplus(:,:,k)=P;
       xplus(:,k)=x;
       measure=1;
    elseif i< sizegps & isintdiv(k-1,ratio)==1 & isnan(d(i)) 
        i=i+1;
    end
    if measure~=1 % no measurement update (needed for smoother)
        Pplus(:,:,k)=P;  % Pplus=Pminus
        xplus(:,k)=x;         % xplus=xminus
    end
    measure=0;   %reset flag
    %Update output variables
    dk(k)=x(1);
    vk(k)=x(2);
end
display(['    ' num2str(i) ' measurement updates in Kalman Filter'])
%display(['    Accelerometer variance: ' num2str(q)])
%display(['    GPS variance: ' num2str(r)])

%SMOOTHING (Rauch-Trung-Striebel algorithm)

clear x
display('Applying RTS smoother...')
x(1,:)=dk;   %make vector of states
x(2,:)=vk;
k=la;   
s(:,k)=x(:,k);   %Initalize smoothed t.series
while k>1
    k=k-1;
    P1=Pplus(:,:,k);     %covariances for gain
    P2=Pminus(:,:,k+1);
    F=P1*A'*inv(P2);   %Smoother gain
    s(:,k)=xplus(:,k)+F*(s(:,k+1)-xminus(:,k+1));   %Update smoothed series
end
clear dk vk
dk=s(1,:); %Update outputs
vk=s(2,:);