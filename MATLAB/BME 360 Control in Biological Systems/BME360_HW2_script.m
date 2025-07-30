% Problem 1
k = 50; % mmol
x = linspace(0, 100, 1000); 
n_val = 1:10;

y = zeros(length(x), length(n_val));
for i = 1:length(n_val)
    n = n_val(i);
    y(:, i) = x.^n ./ (k^n + x.^n);
end

figure;
hold on;
for i = 1:length(n_val)
    plot(x, y(:, i), 'LineWidth', 2);
end
hold off;

title('Hill Function for Different Values of n');
xlabel('Concentration of x (mmol)');
ylabel('Fraction of occupied ligands (y)');
legend('n=1', 'n=2', 'n=3', 'n=4', 'n=5', 'n=6', 'n=7', 'n=8', 'n=9', 'n=10', "Location", "best");
grid on;

saveas(gcf, 'hill_function_plot.png');

% Problem 3
tspan = 1:100
r_val = [-1/2, -1/8, 1, -1/8, -1/2]

r1 = linspace(r_val(1),r_val(1), 20)
r2 = linspace(r_val(2),r_val(2), 20)
r3 = linspace(r_val(3),r_val(3), 20)
r4 = linspace(r_val(4),r_val(4), 20)
r5 = linspace(r_val(5),r_val(5), 20)

r_group = [r1, r2, r3, r4, r5]

[t1, x1] = ode23(@(t,x) r(1)*x + x.^3 - x.^5, 1:20, 1);
[t2, x2] = ode23(@(t,x) r(2)*x + x.^3 - x.^5, 21:40, x1(20,1));
[t3, x3] = ode23(@(t,x) r(3)*x + x.^3 - x.^5, 41:60, x2(20,1));
[t4, x4] = ode23(@(t,x) r(4)*x + x.^3 - x.^5, 61:80, x3(20,1));
[t5, x5] = ode23(@(t,x) r(5)*x + x.^3 - x.^5, 81:100, x4(20,1));

x_group = [x1', x2', x3', x4', x5']

subplot(2,1,1);
plot(tspan, r_group);
title("Hysteresis: r vs. t")
xlabel("t")
ylabel("r")

saveas(gcf, 'hyst_r_vs_t_plot.png');

subplot(2,1,2);
plot(tspan, x_group);
title("Hysteresis: x vs. t")
xlabel("t")
ylabel("x")

saveas(gcf, 'hyst_x_vs_t_plot.png');