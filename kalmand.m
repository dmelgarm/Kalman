function [dk vk]=kalmand(d,a,Td,Ta,q,r)

% D.Melgar 07/2010
% 
% Perform Kalman filter on accelerometer and GPS data to obtain 'true' 
% displacements and velocities (see Smyth & Wu, 2006). The time series have
% to be aligned on the FIRST sample and da/dd must be an integer. The
% acceleration time series can be longer than the gps time series because 
% position will jsut be computed with time update but NOT viceversa.
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
%P=[1 0;0 1];
P=[0 0;0 0];
%Initialize Kalman gain matrix
K=[0 0;0 0];
%Initialize State-Space matrices
A=[1 Ta;0 1];
B=[Ta^2/2;Ta];
H=[1 0];


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
sizegps=size(d,1);
for k=1:la
    %Time update
    Q=[q(k)*Ta^3/3 q(k)*Ta^2/2;q(k)*Ta^2/2 q(k)*Ta]; 
    %predict state
    x=A*x+B*a(k);
    %predict covariance
    P=A*P*A'+Q;
    %Measurement update (if GPS is available)
    if i< sizegps & isintdiv(k-1,ratio)==1  & ~isnan(d(i)) 
        R=r(i)/Td;
        %compute Kalman gain
        K=P*H'/(H*P*H'+R);
        %read in GPS displacement
        %update state
        x=x+K*(d(i)-H*x);
        %update covariance (Joseph stabilized version)
        P=(eye(2)-K*H)*P*(eye(2)-K*H)'+K*R*K';
        i=i+1;
    elseif i< sizegps & isintdiv(k-1,ratio)==1 & isnan(d(i)) 
        i=i+1;
    end
    %Update output variables
    dk(k)=x(1);
    vk(k)=x(2);
end
display(['    ' num2str(i) ' measurement updates in Kalman Filter'])
% display(['    Accelerometer variance: ' num2str(q)])
% display(['    GPS variance: ' num2str(r)])
dk=dk';
vk=vk';
