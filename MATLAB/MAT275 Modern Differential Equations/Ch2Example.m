%Create a function Ch2Example.m and enter the following commands
%function dy= Ch2Example(x,y)
%
%The original ode is dy/dx=1/(x∧2+y∧2)
%
%dy = 1/(x∧2 + y∧2);
%end of function Ch2Example.m
function dy=Ch2Example(x,y)
dy = 1/(x.^2 + y.^2);
end