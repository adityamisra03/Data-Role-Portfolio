
%% Lecture 13


%% Interpolation polynomials

%Taylor error
clear;
NMAX=20;
x=1;x0=0;
testvals=zeros(NMAX,1);
for i=1:NMAX
    testvals(i)=polye(x,x0,i);
end
subplot(1,2,1);plot(testvals)
hold on;
set(gca,'Fontsize',16);
plot(ones(NMAX,1)*exp(1),'r*')
xlabel('Number of Polynomial Terms in Taylor Series');
ylabel('Estimate of e');

subplot(1,2,2);
plot(testvals-ones(NMAX,1)*exp(1),'r*');
set(gca,'Fontsize',16);
xlim([14 16]);
xlabel('Number of Polynomial Terms in Taylor Series');
ylabel('Error in estimate of e');


%% Interpolation first example
clear
A=[1 1 1 1;
    1 2 4 8;
    1 3 9 27;
    1 5 25 125];
RHS=[1.06; 1.12; 1.34; 1.78];

coeffs=A\RHS; %coeffs organized as [d; c; b; a] relative to example
val4=polyval(flip(coeffs),4)
figure;
plot([1 2 3 5],RHS,'b*');
hold on
stem(4,val4);
set(gca,'Fontsize',16);
xlabel('x value')
ylabel('y value')
legend('input values','interpolated value','Location','Southwest')

%% Interpolation second example
x=(0:0.2:1)';
y=log(1+x);
M=zeros(6,6);
M(1,:)=[0 0 0 0 0 1];
for i=5:-1:1
    M(2:6,i)=x(2:6).^i;
end
M(:,6)=1;
coeffs=M\y;

%% ln(1+x) example using Taylor Series
clear;
x=[-1:0.01:1];

T=Tln(x);
P=testpoly(x);%testpoly is just a polynomial containing the coefficients from the last example

figure
plot(x,T);hold on;
plot(x,P,'r');
set(gca,'FontSize',16);
legend({'Taylor Series Approximation N=5';'Test Polynomial'},'Location','southeast');

actual=log(1+x);
plot(x,actual,'g')
stem(0,log(1+0));
ylabel('Approximation to log(1+x)');
xlabel('x')

legend({'Taylor Series Approximation N=5';'Test Polynomial';'Actual';'Pivot'},'location','southeast');

%now look at the error
figure;
plot(x,T-actual);
set(gca,'FontSize',16);
hold on;
plot(x,P-actual,'r')

xlim([0 1]);
ylabel('Error in estimations of log(1+x)');
xlabel('x')
legend({'Taylor Series Approximation N=5';'Test Polynomial';},'Location','northwest');


%% Lagrange interpolation polynomial
clear
x=[0 0.4 0.8 1.2];
% 4 points, can define L3, for k=1, 2, 3
y=cos(x);

myx=pi/4;

myy=L3k(myx,x,3)*y(4)+L3k(myx,x,2)*y(3)+L3k(myx,x,1)*y(2)+L3k(myx,x,0)*y(1);

%% x^3 example
%linear
clear;

N=1;
x=[-1 0];
y=x.^3;

myx=-0.5;
myy=myx^3;

myylinear=0;
for i=0:N
    myylinear=myylinear+Lk(myx,x,N,i)*y(i+1);
end

%quadratic
N=2;
x=[-1:1];%this makes three points
y=x.^3;

myyquadratic=0;
for i=0:N
    myyquadratic=myyquadratic+Lk(myx,x,N,i)*y(i+1);
end

%cubic
N=3;
x=[-1:2];
y=x.^3;

myycubic=0;
for i=0:N
    myycubic=myycubic+Lk(myx,x,N,i)*y(i+1);
end

%based on different points
N=2;
x=1:3;
y=x.^3;
myyquadraticalt=0;
for i=0:N
    myyquadraticalt=myyquadraticalt+Lk(myx,x,N,i)*y(i+1);
end

N=2;
x=[0:2];
y=x.^3;
myyquadraticaltb=0;
for i=0:N
    myyquadraticaltb=myyquadraticaltb+Lk(myx,x,N,i)*y(i+1);
end

errorlin=myy-myylinear;
errorquadalt=myy-myyquadraticalt;
errorcubic=myy-myycubic;

%% Lagrange vs Newton Polynomials

%Lagrange Interpolation polynomials are sums of terms in x^N, no lower
%order terms.  Lots of computation in the (x-x)(x-x)...(x-x)/(x-x)(x-x)
%terms.

%tire data: first column is tire width in mm, second column is final wheel
%circumference. 700 mm wheels assumed.
clear

data=[23	2105
25  2105
28	2143
32	2160
38	2184
40 2200
45 2242];
tirewidthvals=data(:,1);
wssvals=data(:,2);
ninputpoints=length(tirewidthvals);

%% Lagrange Interpolation polynomial
tirewidthfine=[15:0.1:49];
napprox=length(tirewidthfine);
lout=zeros(size(tirewidthfine));
for i=1:napprox
    for j=1:ninputpoints
        lout(i)=lout(i)+wssvals(j)*Lk(tirewidthfine(i),tirewidthvals,ninputpoints-1,j-1);
    end
end

plot(tirewidthvals,wssvals,'m*')
hold on;

plot(tirewidthfine,lout);
set(gca,'Fontsize',14)
xlabel('Tire width (mm)')
ylabel('Tire working circumference (mm)')
ylim([2000 2500]);
%% Newton Polynomial

%Evaluate and plot Newton Polynomial for same data as above
X=tirewidthvals;Y=wssvals;
n=length(tirewidthvals);
D(:,1)=Y';
for j=2:n
    for k=j:n
        D(k,j)=(D(k,j-1)-D(k-1,j-1))/(X(k)-X(k-j+1));
    end
end
A=diag(D);


nout=zeros(size(tirewidthfine));
for i=1:length(tirewidthfine)
    nout(i)=A(1,1);
    term=1;
    for j=0:n-2
        term=term*(tirewidthfine(i)-tirewidthvals(j+1));
        nout(i)=nout(i)+A(j+2)*term;
    end
end
        
% compare the two expressions directly
%figure;
plot(tirewidthfine,nout,'k')
set(gca,'FontSize',14);
legend({'Anchor Points';'Lagrange Interpolation Polynomial';'Newton Interpolation Polynomial'},'Location','North');

%They match at the anchor points and overwrite each other

%% Polynomial Wiggle
%% Lagrange interpolation polynomial
%% Polynomial Wiggle
%% Lagrange interpolation polynomial
clear
close all
x=[0:0.2:1.2];%[0:2e-2:1.2];
N=length(x)-1;
y=cos(x);

myx=pi/4;
mysum=0;
for i=0:N
    myy=mysum+Lk(myx,x,N,i)*y(i+1);
end
%plot this
myfinex=[0:1e-3:1.2];
myfiney=zeros(size(myfinex));
for j=1:length(myfinex)
    myfine=myfinex(j);
    mysum=0;
    for i=0:N
        mysum=mysum+Lk(myfine,x,N,i)*y(i+1);
    end
    myfiney(j)=mysum;
end
figure;
plot(x,y,'r*');
set(gca,'FontSize',16)
hold on;
plot(myfinex,myfiney,'b');
xlabel('x')
ylabel('Approximation Value');
legend('Coarse','Fine')