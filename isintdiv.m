function f=isintdiv(a,b)

%DMM 06/2010
%Decide if a/b produces an integer
%Return 1 if true 0 if false

if a/b - floor(a/b) > 0
    f=0;
else
    f=1;
end