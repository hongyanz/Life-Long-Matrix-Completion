%% MATLAB code for life-long matrix completion for bounded noise
% M: underlying matrix to be recovered
% d: sampling parameter
% eps: bounded noise level
function [L_hat, U_hat, basis_index] = mc_bo(M, d, eps)


%% Initialization
[m, n] = size(M);
U_hat = M(:, 1);
Omega = randperm(m)';
Omega = Omega(1:d);
L_hat = zeros(m, n);
L_hat(:, 1) = M(:, 1);
k = 1;
eta = sqrt(d*k*eps/m)/1.6;
basis_index = 1;

%% Start
for t = 2:n
    if norm(M(Omega, t)-U_hat(Omega, :)*pinv(U_hat(Omega, :))*M(Omega, t)) >= eta
        U_hat = [U_hat, M(:, t)];
        basis_index = [basis_index, t];
        L_hat(:, t) = M(:, t);
        Omega = randperm(m)';
        Omega = Omega(1:d);
        k = k+1;
        eta = sqrt(d*k*eps/m)/1.6;
    else
        L_hat(:, t) = U_hat*pinv(U_hat(Omega, :))*M(Omega, t);
    end
end