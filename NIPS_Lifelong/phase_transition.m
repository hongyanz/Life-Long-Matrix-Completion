clc;
clear;

addpath PROPACK;

m = 50;
n = 500;
success = zeros(50, 50);

for times = 1:1
    times
    for r = 1:1:m
        r
        for d = 1:1:m
            d
            M = rand(m, r)*rand(r, n);
            Omega = randperm(m*n);
            Omega = Omega(1:d*n);
            [L_hat, su] = mc_sp(M, d, r);
            %[L_hat, ~, ~] = unobs_RPCA(M, inf, Omega);
%             if norm(L_hat-M) < 1e-5
%                 success(d, r) = success(d, r)+1;
%             end
            if su == 1
                success(d, r) = success(d, r)+1;
            end
        end
    end
end