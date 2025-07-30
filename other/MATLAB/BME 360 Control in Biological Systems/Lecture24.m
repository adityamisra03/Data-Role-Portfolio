%% Lecture 24

%% Ordinary Differential Equations
clear
clc
close all

%% Compound Interest Example
%Euler's method
%y'=My
%Euler's method yields y_(k+1)=y_k+Ry_k h = y_k(1+Rh)
%Assume the Interest yield is to be calculated over 20 years with an
%annual interest rate of 10% The initial investment is $1000
R=0.1;
Td=20*365;
Tm=20*12;
Tyr=20;
%compute stepsizes, in units of years
hd=1/365;
hm=1/12;
hyr=1;

yyr=zeros(Tyr,1);
ym=zeros(Tm,1);
yd=zeros(Td,1);

%First way, calculate in years
yyr(1)=1000;
for i=2:Tyr
    yyr(i)=yyr(i-1)*(1+R*hyr);
end

%Now in months
ym(1)=1000;
for i=2:Tm
    ym(i)=ym(i-1)*(1+R*hm);
end

%Finally, in days
yd(1)=1000;
for i=2:Td
    yd(i)=yd(i-1)*(1+R*hd);
end

%plot
plot(hd:hd:Td/365,yd,'bo-');
hold on;
plot(hm:hm:Tm/12,ym,'r-*');
plot(hyr:Tyr,yyr,'k--+');
set(gca,'Fontsize',14);
legend({'Days';'Months';'Years'},'Location','NorthWest');
xlabel('Time (years)');
ylabel('Investment value ($)');
%% Change range
xlim([5 7]);
xlim([5 5.5])

%% Taylor series method
% see video
%y'=t^2-y with y(0)=1
%tk=0 tk+1=0.2
%yk=1
%yk+1=yk+d1*h+1/2*d2*h^2+1/6*d3*h^3

h=0.2;
d1=0-1;
d2=2*0-1*(0^2-1);
d3=2-2*0+(0^2-1);
yk3=1+d1*h+0.5*d2*h^2+d3*h^3/6 %N=3
yk2=1+d1*h+0.5*d2*h^2 %N=2

%compare with actual value at t=0.2
ykactual=-exp(-0.2)+0.2^2-2*0.2+2
fprintf('Difference between N=3 Taylor and actual is %f\n',yk3-ykactual);
fprintf('Difference between N=2 Taylor and actual is %f\n',yk2-ykactual);

%% Taylor series method for y'=t^2-y and y(0)=1, summary
%yk+1=yk+(t^2-yk)h+(2t-(t^2-yk)h^2/2+2h^3/6
h=0.2;
timepoints=[0:0.2:1];
y=zeros(size(timepoints));
y(1)=1;
for k=1:5
    t=timepoints(k);
    y(k+1)=y(k)+(t^2-y(k))*h+(2*t-(t^2-y(k)))*h^2/2;%+(2-(2*t-(t^2-y(k))))*h^3/6;
end

%actual solution
actual=-exp(-timepoints)+timepoints.^2-2*timepoints+2;
    
figure;
set(gca,'Fontsize',14);
subplot(1,2,1);
plot(timepoints,y)
hold on;plot(timepoints,actual,'r')
xlabel('Time');
ylabel('y (N=2) or exact solution');
legend({'y (N=2)';'exact solution'});

subplot(1,2,2);plot(timepoints,y-actual);
xlabel('Time');
ylabel('Difference between y (N=2) and actual');


%% Predictor Corrector methods
%Heun method for y'=-y^2
yHeun=zeros(4,1);
h=0.1;
yHeun(1)=1;%value at t=0

%now guess value at t=0.1
%Euler predictor
yHeun(2)=yHeun(1)-h*yHeun(1)^2;
%Corrector
yHeun(2)=yHeun(1)-h/2*(yHeun(1).^2+yHeun(2).^2);

%move to next point t=0.2

%Euler predictor
yHeun(3)=yHeun(2)-h*yHeun(2)^2;
%Corrector
yHeun(3)=yHeun(2)-h/2*(yHeun(2).^2+yHeun(3).^2);


%make a loop
nmax=21;
yHeun(1)=1;
h=0.1
t(1)=0;
for n=2:nmax
    t(n)=t(n-1)+h;%update time 
    yHeun(n)=yHeun(n-1)-h*yHeun(n-1).^2;
    yHeun(n)=yHeun(n-1)-h/2*(yHeun(n-1).^2+yHeun(n).^2);
end


%Another example of predictor corrector method to solve y'=-y^2
%clear
yPC=zeros(nmax,1);
h=0.1;
yPC(1)=1;
yPC(2)=0.909091;
%predictor
yPC(3)=yPC(1)-2*h*yPC(2)^2 %weights are 1 and -2
%corrector
yPC(3)=yPC(2)-h/2*(yPC(3)^2+yPC(2)^2);
%try the corrector again
yPC(3)=yPC(2)-h/2*(yPC(3)^2+yPC(2)^2);

%Run forward to next value
%predictor
yPC(4)=yPC(2)-2*h*yPC(3)^2
%corrector
yPC(4)=yPC(3)-h/2*(yPC(4)^2+yPC(3)^2);
%run it again
yPC(4)=yPC(3)-h/2*(yPC(4)^2+yPC(3)^2);


%make it a loop
yPC(1)=1;
yPC(2)=0.909091;
for n=3:nmax
    t(n)=t(n-1)+h;
    yPC(n)=yPC(n-2)-2*h*yPC(n-1)^2;
    yPC(n)=yPC(n-1)-h/2*(yPC(n)^2+yPC(n-1)^2);
    yPC(n)=yPC(n-1)-h/2*(yPC(n)^2+yPC(n-1)^2);
end
figure
plot(t,yHeun,'-*')
set(gca,'FontSize',14);
hold on;
plot(t,yPC,'r')
xlabel('Time');
ylabel('Solution to ODE, y''=-y^2');
legend('Heun','Alternative PC method');
hold on;
tExact=t;
plot(tExact,1./(tExact+1),'g');
legend('Heun','Alternative PC method','Exact');

figure
plot(t,yHeun'-1./(tExact+1),'b');
set(gca,'FontSize',16)

hold on;
plot(t,yPC'-1./(tExact+1),'r');
legend('Heun','Alternative PC method');
ylabel('Error')

%% Enzyme Kinetics
clear
p.k1=1e-6;
p.km1=1e-8;
p.k2=1e-2;
p.tf=1000;
p.S0=100;
p.E0=10;
p.ES0=0;
p.P0=0;
y0=[p.S0 p.E0 p.ES0 p.P0];

options = odeset('AbsTol', 1e-9, 'RelTol', 1e-6);
[t y] = ode15s(@enzymerhs, [0 p.tf], y0, options, p);

figure;
subplot(2,2,1);
plot (t, y(:,1)); 
set(gca,'FontSize',14);
xlabel ('Time (s)'); ylabel ('S'); title ('Concentration profile ');
subplot(2,2,2);
plot (t, y(:,2));
set(gca,'FontSize',14);
xlabel ('Time (s)'); ylabel ('E'); title ('Concentration profile');
subplot(2,2,3);
plot (t, y(:,3)); 
set(gca,'FontSize',14);
xlabel ('Time (s)'); ylabel ('ES'); title ('Concentration profile');
subplot(2,2,4);
plot (t, y(:,4)); 
set(gca,'FontSize',14);
xlabel ('Time (s)'); ylabel ('P'); title ('Concentration profile');


%% Runge Kutta Pendulum
%Based on EulerPendulum by Kevin Berwick
clear
close all
length= 1;
%pendulum length in metres
g=9.8;
% acceleration due to gravity
npoints = 2500;
%Discretize time into 250 intervals
dt = 0.04;
% time step in seconds
omega = zeros(npoints,1);
% initializes omega, a vector of dimension npoints X 1,to being all zeros
theta = zeros(npoints,1);
% initializes theta, a vector of dimension npoints X 1,to being all zeros
time = zeros(npoints,1);
% this initializes the vector time to being all zeros
theta(1)=0.2;
% you need to have some initial displacement, otherwise the pendulum will
%not swing

for step = 1:npoints-1
% loop over the timesteps using Euler Method
theta(step+1)=theta(step)+omega(step)*dt;
omega(step+1)=omega(step)-g/length*theta(step)*dt;
time(step+1)=time(step)+dt;
end
plot(time,theta,'r');
set(gca,'FontSize',16);
%plots the numerical solution in red
xlabel('time (seconds) ');
ylabel('theta (radians)');
ylim([-pi pi]);
xlim([0 10]);

%another method ('Euler Cromer Method');
omega_EC=zeros(npoints,1);
theta_EC=zeros(npoints,1);
theta_EC(1)=0.2;
for step = 1:npoints-1
% loop over the timesteps
theta_EC(step+1)=theta_EC(step)+omega_EC(step)*dt;
omega_EC(step+1)=omega_EC(step)-g/length*theta_EC(step+1)*dt;
time(step+1)=time(step)+dt;
end
hold on;
plot(time,theta_EC,'k');
%plots the numerical solution in red
xlabel('time (seconds) ');
ylabel('theta (radians)');
ylim([-pi pi]);
xlim([0 10]);
legend({'Euler';'EC';})


%use RK method N=2 (Heun)
omega_rk2=zeros(npoints,1);
theta_rk2=zeros(npoints,1);
theta_rk2(1)=0.2;
for step = 1:npoints-1
% loop over the timesteps
theta_dash_rk2=theta_rk2(step)+0.5*omega_rk2(step)*dt;
omega_dash_rk2=omega_rk2(step)-0.5*theta_rk2(step)*dt;
theta_rk2(step+1)=theta_rk2(step)+omega_dash_rk2*dt;
omega_rk2(step+1)=omega_rk2(step)-g/length*theta_dash_rk2*dt;
time(step+1)=time(step)+dt;
end
plot(time,theta_rk2,'g');
%plots the numerical solution in red
xlabel('time (seconds) ');
ylabel('theta (radians)');
ylim([-pi pi]);
xlim([0 10]);
legend({'Euler';'EC';'RK2';})

%Implement solution using ode45 (RK45);
%define function pendulum
[time_RK4,thetaomega_RK4]=ode45(@pendulum,[0 10],[0.2;0]);
plot(time_RK4,thetaomega_RK4(:,1),'m');
xlabel('time (seconds) ');
ylabel('theta (radians)');
ylim([-pi pi]);
xlim([0 10]);
legend({'Euler';'EC';'RK2';'RK4';})