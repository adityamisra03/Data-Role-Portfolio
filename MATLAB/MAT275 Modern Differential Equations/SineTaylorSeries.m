function f=SineTaylorSeries(N,x1)
%This function takes two inputs N and x1. It compares the Taylor
%expansion of sine with (N+1) terms to the actual sine over
%the interval [-x1,x1].
x=-x1:0.1:x1;
SinePlot=zeros(size(x));
for k=0:N
SinePlot = SinePlot + (-1)^k*x.^(2*k+1)/factorial(2*k+1);
end
plot(x,SinePlot,'b')
hold on
plot(x,sin(x),'r+')
hold off
end