%% Assemble:
%Find: the mass flow rates and velocities of flow at the end of the arteriole bifurcation.
 
%Diagram:
%%
% 
% <<C:\Users\Public\Pictures\problem 3 1.PNG>>
% 
 
%% Assemble:
%   Assume:  Open, Non-reactive, steady state
%   Basis: 
 
 
    v = [20 18];                    % in cm/s
    
% Given data
 
 
    D = [0.2 0.17 0.15 0.12];      % in cm
    rho = 1.056;                    % in gm/cm^3
 
%% Calculate
 
    A = pi * (D/2).^2;
    m(1) = v(1)*A(1)*rho;
    m(2) = v(2)*A(2)*rho;
    
    M = [ 1 0 0; 0 1 0; 1 -1 -1];
    b = [m(1); m(2); 0];

    x = M\b;
    x(4) = x(2)/2;
    x(6) = x(3)/2;
    
    v(4) = x(4)/(A(3)*rho);
    
    v(6) = x(6)/(A(4)*rho);
    
    
    %% Finalize
    % 
    %  
    
    fprintf('the mass flow rate and velocity in branches 4 and 5 are \n')
    fprintf('%g gm/s and %f cm/s, respectively \n\n', x(4), v(4))
    
    fprintf('the mass flow rate and velocity in branches 6 and 7 are \n')
    fprintf('%g gm/s and %f cm/s, respectively', x(6), v(6))
    
    
    
    
    