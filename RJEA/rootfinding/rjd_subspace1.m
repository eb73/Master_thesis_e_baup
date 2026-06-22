function our_solutions = rjd_subspace1(A, nb_LC)
    d = size(A,1);
    n = size(A{1},1);     % size of commuting matrices 
    
    Xs = cell(nb_LC, 1);
    Ys = cell(nb_LC, 1);
    for iter=1:nb_LC
        mu = randn(d,1) + 1i*randn(d,1);
        mu = mu/norm(mu);
        M =  mu(1) * A{1};
        for j = 2:d
            M = M + mu(j)*A{j};
        end
        [X,~,Y] = eig(M);
        Xs{iter} = X;
        Ys{iter} = Y;
    end

    matching = greedyMatching(Xs);
    for t=2:nb_LC
        Xs{t} = Xs{t}(:, matching{t-1});
    end

    X = zeros(n, n);
    %Y = zeros(n, n);
    for idx=1:n % iterates on the pairs X1(:, i) - X2(:, i)
        % extract eigenvectors
        U = zeros(n, nb_LC);
        for col_u = 1:nb_LC
            U(:, col_u) = Xs{col_u}(:, idx);
        end
        [U, ~] = qr(U,"econ"); % U = orth(U); 
        B = cell(d,1);
        for j = 1:d
            B{j} = U' * A{j} * U; % nb_LCxnb_LC matrices
        end

        nB = size(B{1},1); % this should be nb_LC
        warning on
        if nB~=nb_LC
            warning(['size of family B should be ', num2str(nb_LC), '*', num2str(nb_LC), ', constructed: ', num2str(size(B{1},1)), '*', num2str(size(B{1}, 2))])
        end
        % search best joint eigenvector in C^nb_LC
        [best_z, ~] = RJEAbestVec1(B, true, 10);
        X(:, idx) = U*best_z;
        X(:, idx) = X(:, idx) / norm(X(:, idx));
    end

    D = cell(d,1);
    D{1} = diag(X'*(A{1}*X));
    our_solutions = zeros(size(A{1},1),d-1);
    for k = 2:d
        our_solutions(:,k-1) = diag(X'*(A{k}*X));
        D{k} = our_solutions(:,k-1);
    end
    for k = 1:d-1
        our_solutions(:,k) = D{k+1};
    end
end