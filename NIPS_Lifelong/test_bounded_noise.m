clc;
clear;
m = 100;
n = 2000;
r = 5;
d = min(ceil(100*r*log(r)), m);
error = zeros(1, n);

%L = rand(m, r)*rand(r, n);
basis = randn(m, r);
L = [basis(:, 1)*ones(1, 200), basis(:, 2), (basis(:, 1)+basis(:, 2))*ones(1, 199), basis(:, 3), (basis(:, 1)+basis(:, 2)+basis(:, 3))*ones(1, 199), basis(:, 4), (basis(:, 1)+basis(:, 2)+basis(:, 3)+basis(:, 4))*ones(1, 199), basis(:, 5), (basis(:, 1)+basis(:, 2)+basis(:, 3)+basis(:, 4)+basis(:, 5))*ones(1, 1199)];
L = normc(L);
E = 0.05*randn(m, n);
M = L+E;
eps = 0;
for i = 1:n
    if eps < norm(E(:, i));
        eps = norm(E(:, i));
    end
end

[L_hat, U_hat, basis_index] = mc_bo(M, d, eps);
rank(L_hat)
basis_index

for i = 1:n
    error(i) = norm(L_hat(:, i)-L(:, i));
end
plot(1:n, error, '.');
hold on;
pre_error = 0.59*[m*sqrt(1*eps)*ones(1, 200)/d, m*sqrt(2*eps)*ones(1, 200)/d, m*sqrt(3*eps)*ones(1, 200)/d, m*sqrt(4*eps)*ones(1, 200)/d, m*sqrt(5*eps)*ones(1, 1200)/d];
plot(1:n, pre_error, 'r-.', 'LineWidth', 1.2);
axis([0, 2000, 0, 3]);
h = legend('Real Error', 'Estimated Error');
set(h, 'Fontsize', 13);
xlabel('Column Index', 'fontsize', 15);
ylabel('Error', 'fontsize', 15);
