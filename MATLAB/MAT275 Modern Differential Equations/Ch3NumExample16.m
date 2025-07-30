function dy= Ch3NumExample16(x,y)
%
%The original ode is y''+y'+y=0

dy=zeros(2,1); %dy is a column vector!
dy(1) = y(2);
dy(2) = -y(2)-y(1);