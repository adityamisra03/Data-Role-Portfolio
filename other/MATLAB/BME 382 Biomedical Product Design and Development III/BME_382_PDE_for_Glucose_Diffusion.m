%% Solving a PDE for Glucose Diffusion
% code by Jade Lariviere | BME 382 S2024
%% READ ME!
% This code uses the MATLAB pdepe function to solve for the concentration
% of glucose in an unsteady state. The equation takes a specific form in
% the script, and any PDE you attempt to solve must have the proper form.
% pdepe function call is: sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan),
% where:
   % m = symmetry constant (cartesion, cylindrical, spherical)
   % pdefun = definition of the equation to solve
   % icfun =  initial conditions
   % bcfun = boundary conditions
   % xmesh = vector of spatial values
   % tspan = vector of time values
% For more detailed documentation on pdepe with examples, type 'help pdepe'
% into your command window or highlight function and press F1!
%% Solver
clear; close all;
%%%%%%%%%%%% fill in your name and student ID as strings here. %%%%%%%%%%%%
name = 'Aditya Misra';
ID = '1221429625';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------- determine units you'd like to use (write as string) --------- %
u_t = 'seconds'; % time
u_c = 'mmol/L'; % concentration, eg. mmol/L, mg/dL
u_d = 'mm'; % distance
% ----------------------------------------------------------------------- %
% let's fill out some basic values controlling the glucose PDE
% helpful conversions:
   % 1 mmol/L = 18.0182 mg/dL
   % 1 mol glucose = 180.156 g
m = 1; % symmetry
L = 1; % distance
t = 100; % time
tspan = linspace(0,t,100); % choose resolution over time
xspan = linspace(0,L,100); % choose resolution over space
%%%% Go add in the rest of the terms in the following subfunctions! :) %%%%
%%%%%%%%% The hashtags should be replaced with your known values. %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Make sure your units cancel! %%%%%%%%%%%%%%%%%%%%%%%
sol = pdepe(m, @pdeGluc, @icGluc, @bcGluc, xspan, tspan); % fill out to solve PDE
% ---------------------- let's plot your solution! ---------------------- %
figure(1); surf(tspan,xspan,sol);
xlabel('Distance (mm)');
ylabel('Time (seconds)');
zlabel('C_{glucose} (mmol/L)');
title('Unsteady State Glucose Concentration Profile');
subtitle(name + " | " + ID);
% let's create a color map of your 1D concentration gradient!
figure(2); surf(xspan,tspan,sol);
xlabel('Distance (mm)'); ylabel('Time (seconds)');
title('Unsteady State Glucose Concentration: 2D');
subtitle(name + " | " + ID);
view(2); colorbar;
%% Our Glucose Subfunctions
% be sure all your units match!
% defined PDE -------------------------------------------------------------
function [c,f,s] = pdeGluc(x,t,u,dudx)
   Rgluc = 0.1; % reaction rate, concentration/s
   Dgluc = 6.7*10^-6; % diffusivity, distance^2/s
   c = 1; % unsteady-state
   f = Dgluc*dudx; % diffusivity * second partial deriv
   s = -Rgluc*u; % source term, assumed negative; ∝ glucose concentration
end
% -------------------------------------------------------------------------
% initial conditions ------------------------------------------------------
function u0 = icGluc(x,CgluI)
   u0 = 100; % concentration of glucose initially in domain
end
% -------------------------------------------------------------------------
% boundary conditions (left and right) ------------------------------------
function [pL,qL,pR,qR] = bcGluc(xL,uL,xR,uR,t)
   % formatted as Cauchy boundary conditions
   CgluI = 0.1; Cglu1= 0; % concentration at boundaries
   pL = uL - CgluI; qL = 0; % terms for left BC, no flux
   pR = uR - Cglu1; qR = 0; % terms for right BC, no flux
end
% -------------------------------------------------------------------------
