%% Basis
strokevolume = 70;  %mL/beat 
heartrate = 60; %beats/min 
rho = 1.056; %gm/mL 

Vdotj = 296; %2AB
Vdoti = Vdotj*rho/rho; %assuming that the density of blood is the same going in and out of the system
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






