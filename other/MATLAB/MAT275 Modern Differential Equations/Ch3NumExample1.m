function dy= Ch3NumExample1(x,y)
%
%The original ode is 2*x*y”+x^2*y’+3*x^3*y=0
%The system of first-order equations is
%u1'=u2
%u2'=(-x/2)*u2-((3*x^2)/2)*u1
%
%We let y(1)=u1, y(2)=u2
%
dy=zeros(2,1); %dy is a column vector!
dy(1) = y(2);
dy(2) = (-1/2)*x*y(2)-(3/2)*x^2*y(1);
%end of function Ch3NumExample1.m
