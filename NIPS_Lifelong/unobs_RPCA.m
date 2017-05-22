function [L_hat, support, iter] = unobs_RPCA(M, lambda, obs)

[m, n] = size(M);
un_obs = setdiff(1:m*n,obs);

rho = 1.05;

Y1 = zeros(m, n);
Y2 = zeros(m, n);
Y3 = zeros(m, n);
J1 = M;
J2 = zeros(m, n);
mu = 0.001;

max_iter = 1000;
max_mu = mu*1e9;
S_hat = zeros(m, n);
L_old = M;
norm_M = norm(M, 'fro');

iter = 0;
tol = 1e-8;
tol2 = 1e-6;
converged1 = false;
converged2 = false;

while ~converged1 || ~converged2
    iter = iter + 1;
    
    %% Update L_hat
    [U, Sigma, V] = svd(J1-Y2/mu);
    for i = 1:min(size(Sigma))
        if Sigma(i, i) > 1/mu
            Sigma(i, i) = Sigma(i, i)-1/mu;
        elseif Sigma(i, i) < -1/mu
            Sigma(i, i) = Sigma(i, i)+1/mu;
        else
            Sigma(i, i) = 0;
        end
    end
    L_hat = U*Sigma*V';
    
    %% Update S_hat
    temp = J2-Y3/mu;
    for i = 1:n
        if norm(temp(:, i))>lambda/mu
            S_hat(:, i) = (norm(temp(:, i))-lambda/mu)/norm(temp(:, i))*temp(:, i);
        else
            S_hat(:, i) = 0;
        end
    end
    
    %% Update J1
    J1(obs) = (Y1(obs)+Y2(obs)-mu*J2(obs)+mu*M(obs)+mu*L_hat(obs))/(2*mu);
    J1(un_obs) = L_hat(un_obs)+Y2(un_obs)/mu;
    
    %% Update J2
    J2(obs) = (Y1(obs)+Y3(obs)-mu*J1(obs)+mu*M(obs)+mu*S_hat(obs))/(2*mu);
    J2(un_obs) = S_hat(un_obs)+Y3(un_obs)/mu;
    
    %% Update Y1, Y2, Y3
    Y1(obs) = Y1(obs)+mu*(M(obs)-J1(obs)-J2(obs));
    Y2 = Y2+mu*(L_hat-J1);
    Y3 = Y3+mu*(S_hat-J2);

    %% one stop criterion is mu_k*||L_old-L_hat||/||M|| < tol2
    converged2 = false;
    stopCriterion2 = mu*norm(L_old - L_hat, 'fro') / norm_M;
    if stopCriterion2 < tol2
        converged2 = true;
    end
    
    %% the other stop criterion is ||P_obs(M-L_hat-S_hat)||/||M|| < tol  
    converged1 = false;
    temp = M - L_hat - S_hat;
    temp(un_obs) = 0;
    stopCriterion1 = norm(temp, 'fro') / norm_M;
    if stopCriterion1 < tol
        converged1 = true;
    end   
    
    mu = min(rho*mu, max_mu);

    L_old = L_hat;

    if (~converged1 || ~converged2) && iter >= max_iter
        disp('Maximum iterations reached') ;
        converged1 = 1;       
        converged2 = 1;       
    end
end

temp = sqrt(sum(S_hat.^2));
support = temp>1e-3*max(temp);