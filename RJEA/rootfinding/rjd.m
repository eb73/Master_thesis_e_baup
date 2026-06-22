function our_solutions = rjd(A)
    d = size(A,1);
    mu = randn(d,1) + 1i*randn(d,1);
    mu = mu/norm(mu);

    M =  mu(1) * A{1};
    for j = 2:d
        M = M + mu(j)*A{j};
    end
    [X,~] = eig(M);
    D = cell(d,1);
    D{1} = diag(X \(A{1}*X));
    our_solutions = zeros(size(A{1},1),d-1);
    for k = 2:d
        our_solutions(:,k-1) = diag(X \(A{k}*X));
        D{k} = our_solutions(:,k-1);
    end
    for k = 1:d-1
        our_solutions(:,k) = D{k+1};
    end
end