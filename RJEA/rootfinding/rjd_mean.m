function [our_solutions,our_solutions2] = rjd_mean(A, nb_LC) 
    d = size(A,1);
    n = size(A{1},1);     % size of commuting matrices 
    tol_zero = 1e-8;
    
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
        Ys{t} = Ys{t}(:, matching{t-1});
    end
    X = zeros(n);
    Y = zeros(n);
    for idx = 1:n
        z1 = Xs{1}(:, idx);
        if abs(norm(z1) - 1) > tol_zero 
            warning("Some column of X1 does not have norm 1")
            z1 = z1/norm(z1);
        end
        x = z1;
        y = Ys{1}(:, idx);
        for t = 2:nb_LC
            zt = Xs{t}(:, idx);
            wt = Ys{t}(:, idx);
            if abs(norm(zt) - 1) > tol_zero
                warning("Some column of Xt does not have norm 1")
                zt = zt/norm(zt);
            end
            % compute angle between z1 & zt
            alpha = z1'*zt;
            % if angle is non zero, align phase of complex vectors
            if abs(alpha) > tol_zero
                mult = exp(-1i*angle(alpha));
                zt = zt*mult; % this is equivalent to z2*conj(alpha)/abs(alpha)
                wt = wt*mult;
            end
            x = x + zt;
            y = y + wt;
        end
        x = x/nb_LC;
        y = y/nb_LC;
        % x is already normalized by contruction, but check in case
        if abs(norm(x) - 1) > tol_zero 
            warning("Needed to normalize constructed eigenvector.")
            x = x/norm(x);
        end
        X(:, idx) = x; 
        Y(:, idx) = y;
    end
    D = cell(d,1);
    D{1} = diag(X'*A{1}*X);
    our_solutions = zeros(size(A{1},1),d-1);
    for k = 2:d
        our_solutions(:,k-1) = diag(X'*A{k}*X);
        D{k} = our_solutions(:,k-1);
    end
    for k = 1:d-1
        our_solutions(:,k) = D{k+1};
    end
    D = cell(d,1);
    D{1} = diag(X\(A{1}*X));
    our_solutions2 = zeros(size(A{1},1),d-1);
    for k = 2:d
        our_solutions2(:,k-1) = diag(X\(A{k}*X));
        D{k} = our_solutions2(:,k-1);
    end
    for k = 1:d-1
        our_solutions2(:,k) = D{k+1};
    end
end