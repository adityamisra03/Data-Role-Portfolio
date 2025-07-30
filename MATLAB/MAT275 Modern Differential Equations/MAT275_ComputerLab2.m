%% Computer Lab 2 Chapter 2 #: 1, 2, 3, 5b, 5d, 6b, 6d, 9
%The "Example 1", "Example 2", and "Example 3" referred to in #1,2,3 are on pp. 115-118 of the textbook.  
% On #6, you may use RK4 or ode45 and choose three initial conditions of your choice.  On #9, print out a table 
% (as shown in the class lecture on 2.5 Numerical Methods) that compares the outputs; 
% you must use Euler.m and RK4.m that you wrote for the previous problems.

%% Chapter 2 problems
fprintf('Chapter 2 \n')

%% 1
fprintf('Problem 1 \n')
clear all;
[X,Y]=meshgrid(-3:.55:3,-3:.55:3); %Note the CAPITAL letters
DY=1./(X.^2+Y.^2); %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2); %the length of each vector
quiver(X,Y,DX./DW,DY./DW,.5,'.'); %plots normalized vectors
xlabel('x');
ylabel('y');

%% 2
fprintf('Problem 2 \n')
x0=-3; xf=3; y0=-3;
options=odeset('refine',10,'AbsTol',1e-15);
[x,y]=ode45(@Ch2Example,[x0,xf],y0,options);
x0=-3; xf=3; y0=-2;
[x1,y1]=ode45(@Ch2Example,[x0,xf],y0,options);
x0=-3; xf=3; y0=-1;
[x2,y2]=ode45(@Ch2Example,[x0,xf],y0,options);
x0=-3; xf=3; y0=0;
[x3,y3]=ode45(@Ch2Example,[x0,xf],y0,options);
hold on %keeps direction field plot open
plot(x,y) %superimpose solutions
axis([-3 3 -3 3])
plot(x1,y1)
plot(x2,y2)
plot(x3,y3)
hold off

%% 3
fprintf('Problem 3 \n')
clear all
close all
x0=0; xf=3.1; y0=1; h=.1;
[x,y]=Euler(@Ch2NumExample,[x0, xf],y0,h);
[x y]
[x1,y1]=RK4(@Ch2NumExample,[x0 xf],y0,h);
[x1 y1]
x2=x0:h:xf;
y2= exp(x2.^2/2); %this is the analytical solution
[x2' y2']
plot(x,y,'b:') %Now graphically compare Euler with RK4
hold on
plot(x1,y1,'k--')
xlabel('x')
ylabel('y')
title('Numerical solns of y\{prime} =xy using EM and RK4 with h=.1, dotted/blue is EM, dashed/black is RK4')
hold off

plot(x1,y1,'k--','LineWidth',2)
hold on
[X,Y]=meshgrid(-3:.3:3,0:.3:4);
DY=X.*Y; %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2);
quiver(X,Y,DX./DW,DY./DW,.4,'.');
xlabel('x');
ylabel('y');
axis([-3 3 0 4])
hold off
%End of Example 3.

%% 5b

fprintf('Problem 5b \n')
clear all;
[X,Y]=meshgrid(-5:.25:5,-5:.25:5); %Note the CAPITAL letters
DY=sin(X.^2); %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2); %the length of each vector
quiver(X,Y,DX./DW,DY./DW,.5,'.'); %plots normalized vectors
xlabel('x');
ylabel('y');

%% 5d
fprintf('Problem 5d \n')
clear all;
[X,Y]=meshgrid(-5:.25:5,-5:.25:5); %Note the CAPITAL letters
DY=1-Y.^2; %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2); %the length of each vector
quiver(X,Y,DX./DW,DY./DW,.5,'.'); %plots normalized vectors
xlabel('x');
ylabel('y');

%% 6b

fprintf('Problem 6b \n')
clear all;
[X,Y]=meshgrid(-5:.25:5,-5:.25:5); %Note the CAPITAL letters
DY=sin(X.^2); %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2); %the length of each vector
quiver(X,Y,DX./DW,DY./DW,.5,'.'); %plots normalized vectors
xlabel('x');
ylabel('y');

options=odeset('refine',10,'AbsTol',1e-15);

x0=-5; xf=5; y0=-2;
[x1,y1]=ode45(@Ch5bNumExample,[x0,xf],y0,options);

x0=-5; xf=5; y0=-1;
[x2,y2]=ode45(@Ch5bNumExample,[x0,xf],y0,options);

x0=-5; xf=5; y0=0;
[x3,y3]=ode45(@Ch5bNumExample,[x0,xf],y0,options);

hold on %keeps direction field plot open
 %superimpose solutions
axis([-5 5 -5 5])
plot(x1,y1)
plot(x2,y2)
plot(x3,y3)
hold off

%% 6d

fprintf('Problem 6d \n')
clear all;
[X,Y]=meshgrid(-5:.25:5,-5:.25:5); %Note the CAPITAL letters
DY=1-Y.^2; %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2); %the length of each vector
quiver(X,Y,DX./DW,DY./DW,.5,'.'); %plots normalized vectors
xlabel('x');
ylabel('y');

options=odeset('refine',10,'AbsTol',1e-15);

x0=-2; xf=2; y0=0;
[x1,y1]=ode45(@Ch5dNumExample,[x0,xf],y0,options);

x0=-2; xf=2; y0=-2;
[x2,y2]=ode45(@Ch5dNumExample,[x0,xf],y0,options);

x0=-2; xf=2; y0=3;
[x3,y3]=ode45(@Ch5dNumExample,[x0,xf],y0,options);

hold on %keeps direction field plot open
%superimpose solutions
axis([-5 5 -5 5])
plot(x1,y1)
plot(x2,y2)
plot(x3,y3)
hold off

%% 9
fprintf('Problem 9 \n')
clear all;
close all;
x0=0; xf=5; y0=-1; h=.1;
[x,y]=Euler(@Ch9NumExample,[x0, xf],y0,h);
[x y]
[x1,y1]=RK4(@Ch9NumExample,[x0 xf],y0,h);
[x1 y1]

syms X Y(X)
eqn=diff(Y,X)==Y(X)+X.^2;
cond=Y(0)==-1;
ysol(X)=dsolve(eqn,cond);
ysol(X)

x2=x0:h:xf; 
y2=exp(x)-2*x-x.^2-2 %this is the analytical solution
%[x2' y2']
plot(x,y,'b:') %Now graphically compare Euler with RK4
hold on
plot(x1,y1,'c--')
xlabel('x')
ylabel('y')
title('')
hold off

plot(x1,y1,'k--','LineWidth',2)
hold on
[X,Y]=meshgrid(-5:.25:5,-5:.25:5);
DY=Y+X.^2; %DY is rhs of original equation
DX=ones(size(DY));
DW=sqrt(DX.^2+DY.^2);
quiver(X,Y,DX./DW,DY./DW,.5,'.');
xlabel('x');
ylabel('y');
axis([-5 5 -5 5])
hold off


fprintf('Problem 9 \n')
clear all;
syms x y(x)
format longg
dy=y(x)+x.^2;

eqn=diff(y,x)==y+x.^2;
cond=y(0)==-1;
ysol(x)=dsolve(eqn,cond);
ysol(x)

[x1,y1]=Euler(@Ch9NumExample,[0,5],-1,0.1);
[x1, y1, double(ysol(x1))]

[x2,y2]=RK4(@Ch9NumExample,[0,5],-1,0.1);
[x2, y2, double(ysol(x2))]

xeulerr4kanalyticalans=[0 -1 0 -1;
        0.1 -1.1 -1.10482895833333 -1.10482908192435;
        1 -2.24688329389 -2.28171685211747 -2.28171817154095;
        5 92.029938167665 111.412882112512 111.413159102577]
difference=[xeulerr4kanalyticalans(:,1) abs(xeulerr4kanalyticalans(:,2)-xeulerr4kanalyticalans(:,4)) abs(xeulerr4kanalyticalans(:,3)-xeulerr4kanalyticalans(:,4))]


