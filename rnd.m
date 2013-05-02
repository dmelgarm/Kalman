function rnum=rnd(num)

%DMM 10/2010
%Round down or up to the nearest integer


rnum=num-floor(num);
if rnum<0.5
    rnum=floor(num);
else
    rnum=ceil(num);
end