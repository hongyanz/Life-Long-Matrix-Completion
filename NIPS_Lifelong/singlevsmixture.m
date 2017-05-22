clc;
clear;
r = 20;
M = randn(100, r)*randn(r, 100);
M(1:30, :) = 0;
[m, n] = size(M);
error = zeros(1, 100);

for times = 1:100
    k = 0;
    for d = 0.01*m : 0.01*m : m
        k = k+1;
        [L_hat, success] = mc_sp(M, d, r);
        if success == 0
            error(k) = error(k)+1;
        end
    end
end

plot(0.01*m : 0.01*m : m, error/times, '.')
axis([0, 100, -0.1, 1.1])