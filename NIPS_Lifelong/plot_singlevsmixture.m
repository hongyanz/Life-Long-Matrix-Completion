clc;
clear;

load('singlevsmixture_single.mat');
error = 200-error;
plot(0.01*m : 0.01*m : m, error/times, '.-r', 'LineWidth', 2.5);
axis([0, 100, -0.1, 1.1]);
hold on;
grid on;

load('singlevsmixture_mixture.mat');
error = 200-error;
plot(0.01*m : 0.01*m : m, error/times, '--', 'LineWidth', 2.5);

LEG = legend('Single Subspace', 'Mixture of Subspaces');
xlabel('Parameter d', 'fontsize', 15);
ylabel('Empirical Probability of Success', 'fontsize', 15);
set(LEG, 'FontSize', 15)