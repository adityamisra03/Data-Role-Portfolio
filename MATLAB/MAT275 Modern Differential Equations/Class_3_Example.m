
%% Integrate
%function tointegrate 

%a = 0;  %start point 
%b = 10; %endpoint  
%Q = integral(@myfun,a,b) 
%end

%function y = myfun(x) 
%y = x.^3; %function to integrate.  Note .^ for power not just ^ 
%end



%% Shortcut
a = 0; %starting point of integration  
b = 4; %ending point of integration 
fun = @(x) x.^3;  
Q = integral(fun,a,b); 
fprintf('The integral of x^3 from 0 to 4 is %2.0f. \n', Q);  
fprintf('This was an example of a short cut to integrate \nwith the integral function. \n'); 
% This is the output of the fprintf's
% The integral of x^3 from 0 to 4 is 64. 
% This was an example of a short cut to integrate
% with the integral function.   
%% Basis
strokevolume = 70;  %mL/beat 
heartrate = 60; %beats/min 
rho = 1.056; %gm/mL 

Vdotj = 296;
Vdoti = Vdotj*rho/rho; 
a = 0; 
b = 30;  % 30 minutes
fun = @(t) Vdotj + 0*t; 
Q = integral(fun,a,b); 
% mdoti - mdotj = 0; 
% Vdoti*rho - Vdotj*rho = 0 
% Vdoti - Vdotj = 0 

%% Finalize
fprintf('the volumetric flow rate (in and out) for left side of the heart is %6.2f mL/min \n',Vdoti) 
fprintf('\n the volume in 30 minutes would be %4.2f L. \n', Q/1000); 






