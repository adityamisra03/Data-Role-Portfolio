function dy= Ch3NumExample2(x,y)
%
%The original ode is y^(4)+x^2*y’+y=cos(x)
%The system of first-order equations is
%u1'=u2
%u2'=u3
%u3'=u4
%u4'=-x^2*u2-u1+cos(x)
%
%We let y(1)=u1, y(2)=u2, y(3)=u3, y(4)=u4
%
dy=zeros(4,1); %dy is a column vector!
dy(1) = y(2);
dy(2) = y(3);
dy(3) = y(4);
dy(4) = -x^2*y(2)-y(1)+cos(x);
%end of function Ch3NumExample2.m