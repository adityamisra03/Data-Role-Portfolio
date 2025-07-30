% Appendix A #: 1,2,3,4,6,7,12,13,18 AND Chapter 1 #: 3,4,19,22,26.
%% Appendix A Problems
fprintf('Appendix A\n')
%% 1
fprintf('Problem 1\n')
clear all;
pi+exp(1)+log(4)
format long
pi+exp(3)+log(1/2)
x=-2:.1:pi/2;
eq1=x.^2-1+sin(x);
plot(x,eq1)
plot(x,eq1,'b--','LineWidth',3)

xlabel('x','Fontsize',12), ylabel('y','Fontsize',12)
legend('x^2-1+sin(x)')
title('A basic plot')
axis([-2 pi/2 -1.5 1.5])
f=sqrt(x+2)+log(x+2); %note that log(x) is ln(x) in MATLAB
plot(x,f)
plot(x,eq1,'b-',x,f,'m.')
xlabel('x','Fontsize',12), ylabel('y','Fontsize',12)
title('Superimposed plots')
legend('x^2-1+sin(x)','sqrt(x+2)+log(x+2)')
axis([-2 1 -1.5 2])
x=0:.1:3;
eq2=log(x);
eq3=exp(-x.^2);
plot(x,eq2,'b')
hold on
plot(x,eq3,'k--')
xlabel('x','Fontsize',12), ylabel('y','Fontsize',12)
axis([0 3 -2 2]), title('ln(x) vs exp(-x^2)')
hold off %anything that follows will not be superimposed
orient tall %makes picture take up full page when you print

%% 2
fprintf('Problem 2 \n')
clear all;
syms x;
f(x)=x^4+2*x^3-8*x^2;
f(3);
subs(f(x),x,3);
eq1=solve(f(x),x);
eq1(4);
factor(f(x));
g(x)=1-x^2;
f(x)*g(x);
factor(f(x)*g(x));
f(x)/x;
simplify(f(x)/x);
h(x)=x^2+2*x+4;
solve(h(x),x);
eq2=solve(h(x),x);
eq2(1);
real(eq2(1));
imag(eq2(1));
2*i+eq2(1);
eval(2*i+eq2(1));
eq1=x^4+3*x^2-7*x+6;
solve(eq1);
%Next 2 lines find roots numerically w/o Symbolic Math Toolbox
p=[1 0 3 -7 6]; % Coefficients of polynomial
roots(p);

%% 3
fprintf('Problem 3 \n')
clear all;
SineTaylorSeries(2,3)
SineTaylorSeries(3,6)
SineTaylorSeries(8,10)

%% 4
fprintf('Problem 4 \n')
clear all;
format long
exp(pi)
pi.^exp(1)
fprintf("e^pi = %12.12f is larger than pi^(e^1) = %12.12f", exp(pi), pi.^exp(1))

%% 6
fprintf('Problem 6 \n')
clear all;
x=-3*pi:0.1:3*pi;
y=sin(2*x);
plot(x,y,"r:")

xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
title('Problem 6');
legend('y=sin(2x)')
orient tall

%% 7
fprintf('Problem 7 \n')
clear all;
x=-1:0.1:3;
y= (1-x).^2-2;
plot(x,y, "c-.")

xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
title('Problem 7');
legend('y=(1-x)^2-2')
orient tall

%% 12
fprintf('Problem 12 \n')
clear all;
x= -2*pi:0.1:2*pi
y1= sin(x)
y2= sin(2*x)
y3= cos(x)
y4= cos(2*x)

subplot(2,2,1),plot(x,y1)
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=sin(x)')

subplot(2,2,2),plot(x,y2)
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=sin(2x)')

subplot(2,2,3),plot(x,y3)
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=cos(x)')

subplot(2,2,4),plot(x,y4)
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=cos(2x)')

%% 13
fprintf('Problem 13 \n')
clear all;
x= 0.1:0.1:10
y= atan(x)
plot(x,y, 'r-')
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=arctan(x)')

hold on
x= 0.1:0.1:10
y= sin(x.^(-1))
plot(x,y, 'c-')
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);
legend('y=arctan(x)','y=sin(1/x)')
hold off

%% 18
fprintf('Problem 18 \n')
clear all;
InfSeries(100);
InfSeries(1000);
InfSeries(10);
InfSeries(000);

%% Chapter 1 Problems

%% 3
fprintf('Problem 3 \n')
clear all;
[X,Y]=meshgrid(-2:.01:2,-2:.01:2);
Z=-1./X+2./X.^2+1./Y-1./Y.^3;
contour(X,Y,Z,[-5 -5])
[C,h]=contour(X,Y,Z,[-5 -5]);
clabel(C,h)
axis([-2 2 0 0.6])
xlabel('x'), ylabel('y')
title('Implicit Plot with C=-5')
[C,h]=contour(X,Y,Z,[-10,2,6,50]);
clabel(C,h)
xlabel('x'), ylabel('y')
title('Implicit Plot with C=-10,2,6,50')

%% 4
fprintf('Problem 4 \n')
clear all;
[x,y]=ode45(@Example4,[0,4],10);
plot(x,y)
subplot(2,1,1),plot(x,y)
title('Numerical Approximation of Soln of (x.^2+1)y{\prime}+4xy=x, y(0)=10')
xlabel('x'), ylabel('y')
x1=0:.1:4;
y1=((1/4)*x1.^4+(1/2)*x1.^2+10)./(x1.^2+1).^2;
subplot(2,1,2),plot(x1,y1,'k')
title('Closed Form Soln of (x.^2+1)y{\prime}+4xy=x,y(0)=10')
xlabel('x'), ylabel('y')

%% 19
fprintf('Problem 19 \n')
clear all;
syms x;
y(x)=sin(x)+2*cos(x);
dy(x)=diff(y(x));
d2y(x)=diff(dy(x));
solution = d2y(x) + y(x);
fprintf('Since the second dervative of y is %s and the y is %s,\n the solution of %s is equal to 0 because they cancel out, so there is a solution \n on (-inf, inf) to d^2y/dx^2 + y = 0. \n', d2y(x), y(x), solution)

%% 22
fprintf('Problem 22 \n')
clear all;
syms x;
y(x)= int(sin(x.^2),x);
x=0:0.01:4;
plot(x,y(x))
xlabel('x','Fontsize',12); 
ylabel('y','Fontsize',12);

%% 26
fprintf('Problem 26 \n')
clear all;
%syms x y(x) C ;

%eq2=y(x)+x*diff(y(x),x)-(y(x)/(2*y(x)-1));
%sol2=(x+(x.^2-C).^(1./2))./(2*x);

clear all;
syms X Y Z;
[X,Y]=meshgrid(-10:.1:10,-10:0.1:10);
Z=(X.^2)-(2.*X.*Y-X.^2).^2;

[C,h]=contour(X,Y,Z,[1,4,8]);
clabel(C,h)
xlabel('x'), ylabel('y')
title('Implicit Plot with C=1,4,8')



