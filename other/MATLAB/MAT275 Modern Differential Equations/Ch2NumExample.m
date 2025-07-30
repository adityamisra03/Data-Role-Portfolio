%Create a function Ch2NumExample.m and enter the following commands
%function dy= Ch2NumExample(x,y)
%
%The original ode is dy/dx=x*y
%
function dy= Ch2NumExample(x,y)
%dy = x*y;
dy=3*x*exp(-y)
%end of function Ch2NumExample.m
end