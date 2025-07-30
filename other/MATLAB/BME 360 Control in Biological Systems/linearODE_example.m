function dxdt = linearODE_example(t,x)
    % This is the place where your function lives and gets called. 
    
    % In order for this function to represent y' = 10*t + 2, see the following:
    
    % dydt corresponds to y' and your desired output. That's why we assign the variable to our ODE.
    % In the case that you want to return more outputs, then you would
    % replace dydt with [something, anotherthing, andwhatever, +...]
    
    % The linearODE_example slot is the name of your function
    
    % The (t,y) is the input arguments. You need both in this case as you have
    % an input for the time interval and for the initial condition of y0
    
    dxdt = x.^2-1;
    %dxdt = sin(x);
end
