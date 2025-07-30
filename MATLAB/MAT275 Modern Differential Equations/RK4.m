%Create a function RK4.m and enter the following commands
%function [xout,yout]=RK4(fname,xvals,y0,h)
%This code implements the 4th order Runge--Kutta method for
%numerically solving the ODE y’=f(x,y)
%
%fname = the function f for the equation.
%x0 =initial x
%n = number of steps to be taken.
%y0 =initial y
%h =stepsize
function [xout,yout]=RK4(fname,xvals,y0,h)
x0=xvals(1); xf=xvals(2);
x=x0; y=y0;
steps=round((xf-x0)/h);
f=feval(fname,x,y);
xout=zeros(steps, 1);% allocates space for xout
yout=zeros(steps, length(f));% allocates space for yout
y=y'; %The ' is needed to match syntax of ode45 in higher dim
xout(1,1)=x; yout(1,:)=y;
for i=1:steps
    k1 = h*f;
    k2 = h*feval(fname,x+(h/2),y+(k1/2));
    k3 = h*feval(fname,x+(h/2),y+(k2/2));
    k4 = h*feval(fname,x+h,y+k3);
    ynext = y +(k1 + 2*k2 + 2*k3 + k4)/6;
    xnext = x+h;
    f = feval(fname,xnext,ynext);
    xout(i,:)=xnext;
    yout(i,:)=ynext'; %we again need the '
    x=xnext;
    y=ynext;
end
end
%end of function RK4.m