function best_sol = rjd_best_res(A, normsA, nb_LC)
    d = size(A,1);
    best_sol = zeros(size(A{1},1),d-1);
    best_res = Inf;
    for iter = 1:nb_LC
        mu = randn(d,1) + 1i*randn(d,1);
        mu = mu/norm(mu);
    
        M =  mu(1) * A{1};
        for j = 2:d
            M = M + mu(j)*A{j};
        end
        [X,~] = eig(M);
        sol_tmp = zeros(size(A{1},1),d);
        for k = 1:d
            sol_tmp(:,k) = diag(X\(A{k}*X));
        end
        res_tmp = compute_residual(A,sol_tmp,normsA); %compute_residual(A, sol_tmp); % here we need to use the same residual as rjd (not poly residual)
        if res_tmp < best_res
            best_res = res_tmp;
            best_sol = sol_tmp(:, 2:d);
        end
    end
end