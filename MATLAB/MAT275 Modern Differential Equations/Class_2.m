% Class 2 Problems

%% Original Problems
%5*u + 6*v + w = 0

%3*v + 2*w = 3

%10*u + 9*v + 3*w = 8

M = [5 6 1; 0 3 2; 10 9 3];
b = [0;3;8];
x = M\b;
fprintf('x =%4.4f \n', x)

%x1^2 - 2*x2 - 3*x3 = 4

%-2*x1 + 2*x2 = 6

%5*x1 + x2 - 7*x3 = -4

syms x1 x2 x3;
E1 = x1^2 - 2*x2 - 3*x3 - 4;
E2 = -2*x1 + 2*x2 - 6;
E3 = 5*x1 + x2 - 7*x3 + 4;
X = solve(E1, E2, E3);
x1 = double(X.x1);
x2 = double(X.x2);
x3 = double(X.x3);
fprintf('x1 = %4.2f ; x2 = %4.2f ; x3 = %4.2f \n', x1, x2, x3)

%% New problems
W = [8 6 7 5; 4 5 2 3; 2 1 2 3];
p = [8;13;3];
d = W\p;
fprintf('d =%4.4f \n', d)

syms v1 v2 v3;
Eq1 = 4*v1^2 - 2*v2 - v3 - 8;
Eq2 = -3*v1 + 2*v2 + 6*v3 - 16;
Eq3 = 4*v1 + 6*v2 - 3*v3 + 20;
V = solve(Eq1, Eq2, Eq3);
v1 = double(V.v1);
v2 = double(V.v2);
v3 = double(V.v3);
fprintf('v1 = %4.2f ; v2 = %4.2f ; v3 = %4.2f \n', v1, v2, v3)

