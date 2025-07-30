clear; clc;

r_values = linspace(0, 3, 1000);
x0 = 0.1;
num_iter = 1000;
threshold = 1e-6;

bifur_data = zeros(length(r_values), num_iter);

for i = 1:length(r_values)
    r = r_values(i);
    for j = 1:num_iter
        f = @(x) r - x - exp(x);
        x_root = fsolve(f, x0);
        x0 = x_root;
        bifur_data(i, j) = x_root;

        if abs(f(x_root)) < threshold
            break;
        end
    end
end

figure;
plot(r_values, bifur_data, '.', 'MarkerSize', 2);
xlabel('r');
ylabel('Fixed points (x)');
title('Bifurcation Diagram');