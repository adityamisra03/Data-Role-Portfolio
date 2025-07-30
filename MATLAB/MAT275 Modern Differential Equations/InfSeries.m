function f=InfSeries(N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sum1 = 0;
for k = 1:N
    sum1 = sum1 + 1/k;
end
fprintf('The sum is %6.6f.\n',sum1)
end
