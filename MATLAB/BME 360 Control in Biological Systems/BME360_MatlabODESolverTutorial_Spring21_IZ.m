%% BME 360 MATLAB ODE SOLVER TUTORIAL - Irene Zhang 1/19/20
% The purpose of the following is to help guide you on using MATLAB's ODE
% solver for a 1D ODE. Information is synthesized from Mathworks MATLAB
% documentation directory

%% ODE Solver - Basic Background

% ODE (Ordinary differential equation) solvers solve differential equations
% using known initial conditions and boundary conditions

% There are different kinds of ODE solvers, the most common ones used being
% ode45 and ode23. Both of these ode solvers are nonstiff, meaning they are
% meant to solve problems that shouldn't be too difficult to evaluate. A
% equation is stiff if it takes too long to solve using nonstiff solvers or
% it will not solve at all. 

% For the sake of this tutorial, we will be using ode23. This solver works great
% for more complex systems and higher order equations as well (the van der Pol equation is a great example of
% a second order ODE, and also plots to be quite beautiful)

%% ODE Solver Example 

% Let's choose the basic linear (1D) ODE  x' = 2*x - 100 as an example
% equation. The time interval will be [0 20] and the initial condition y0 =
% 0. 

% Though this ODE is very simple, it is good practice to save your ODE in a
% separate m file as a function. See the attached linearODE_example.m 

% The syntax of the basic ode23 solver is as follows:
% [t,x] = ode23(odefun, tspan, x0)

% where [t,x] is your outputs of evaluation points (times evaluated at) and solutions (x) 
% odefun is the name of your ODE function(s)
% tspan is the interval of integration (given)
% x0 is your initial condition

% So, based on this information, the first step is to set up your solver by
% coding the given values:

close all;
tspan = [0 1];
x0 = 100;

% Next, call the ode solver and call the ode function you have written

[t,x] = ode23('linearODE_example',tspan, x0);

% We can also print the ode function in this main file for easier viewing
% afterwards or if published:

type linearODE_example.m

% Now let's plot the solution:

plot(t,x)
title('Sample ODE Plot: dxdt = 2x - 100');
xlabel('t');
ylabel('x');

% For simple ODEs like this, you can also skip the separate function and
% write it as an anonymous function within the ode23 call. It will get you
% the same answers! But in general, it can be easier to understand and read
% your code if you make the function separate. 

%% ODE Solver Help 

% A highly underrated command in MATLAB is doc ___. If you type doc
% "whatever you are looking for or confused by," then you will get all the
% official documentation in the MATLAB library. Great for when you are
% troubleshooting! 

% If you have anymore questions or issues on ODE solver or understanding
% how to use it, please feel free to email me (izhang1@asu.edu) or come see
% me in my office hours (Mondays 3-5 pm). Good luck! 