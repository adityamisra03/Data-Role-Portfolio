%Create a function Euler.m and enter the following commands
%function [xout,yout]=Euler(fname,xvals,y0,h)
%This code implements the Euler Method for numerically
%solving y’=f(x,y)
%
%fname=the function f for the ODE
%xvals = vector that contains the initial x0 and xf
%y0 = initial y
%h = stepsize
function [xout,yout]=Euler(fname,xvals,y0,h)
x0=xvals(1); xf=xvals(2);
x=x0; y=y0;
steps=(xf-x0)/h;
xout=zeros(steps+1, 1);% allocates space for xout
yout=zeros(steps+1, 1);% allocates space for yout
xout(1)=x; yout(1)=y;
for j=1:steps
    f=feval(fname,x,y);
    x=x+h;
    y=y+h*f;
    xout(j+1)=x;
    yout(j+1)=y;
end
end
%end of function Euler.m
