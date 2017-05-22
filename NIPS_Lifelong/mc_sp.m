%% MATLAB code for life-long matrix completion for sparse random noise
% M: underlying matrix to be recovered
% d: sampling parameter
% r: rank of M
function [L_hat, success] = mc_sp(M, d, r)


%% Initialization
[m, n] = size(M);
U_hat = M(:, 1);
Omega = randperm(m)';
Omega = Omega(1:d);
L_hat = zeros(m, n);
L_hat(:, 1) = M(:, 1);
success = 0;
k = 0;

%% Start
for t = 2:n
    if norm(M(Omega,t)-U_hat(Omega, :)*pinv(U_hat(Omega, :))*M(Omega,t)) > 1e-5
        k = k+1;
        U_hat = [U_hat, M(:, t)];
        L_hat(:, t) = M(:, t);
        Omega = randperm(m)';
        Omega = Omega(1:d);
    else
        L_hat(:, t) = U_hat*pinv(U_hat(Omega, :))*M(Omega, t);
    end
end

precision = norm(L_hat-M, 'fro');
if precision < 1e-5 && size(U_hat, 2) == r
    success = 1;
end