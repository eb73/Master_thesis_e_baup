function [rq1, rq2, res, err, cnd, normsX, d2] = newRayleighMean(A,N,exact,clus,complex,noise,Aset,nb_LC)
% Computes the approximate joint eigenvalues with RJEA+mean(k=nb_LC)
%
% Input:
%  - A: a set of d commuting n x n matrices (usually perturbet matrices from Aset)
%  - N: number of runs
%  - exact: n x d matrix with exact joint eigenvalues of A for comparison
%  - clus: index of eigenvalue or vector of indices of a cluster
%  - complex: (1) use complex random coefficients, 0: use real
%  - noise: (0) amount of Gaussian noise added to A(mu) before computing eigenvector matrix X
%  - Aset: (A) original matrices used for computing residuals
%  - nb_LC: (2) number of linear combination used to pick best solution
%  
% Output:
%  - rq1: cell of size d with n X N matrices, where (i,j)-th element of k-th
%        cell is 1-sided RQ approximation for exact(i,k) in run j
%  - rq2: same as rq1 when 2-sided RQ is used
%  - res: vector of size N x 2 with residuals
%         v(k) = norm((A{k}-rho(j,k)*I)*X(:,j))/norm(A{k}), 
%         where rho(j,k) is an eigenvalue approximation (RQ1 & 2)
%  - err: N x 4 matrix with 2-norm of the difference between the computed 
%         joint eigenvalue and the exact eigenvalue (RQ1 & 2) (first
%         for cluster and then for the first index from the cluster)
%  - cnd: vector of size N with condition numbers of eigenvector matrix X
%  - d2: values d(mu), s(mu), t(mu) and gap
%  - normsX : values [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)] from diagonalization

class_t = superiorfloat(A{:});

% default values 
if nargin<5, complex = 1; end 
if nargin<6, noise = 0;   end 
if nargin<7, Aset = A; end
if nargin<8, nb_LC = 2; end

n = size(A{1},1);     % size of commuting matrices 
d = length(A);        % number of commuting matrices in set A 
tol_zero = 1e-8;

% we compute complements if indices in clusters A and B
clus_comp = 1:n;
clus_comp(clus) = 0;
clus_comp = find(clus_comp~=0);
clus_comp_size = length(clus_comp);
clus_size = length(clus);
res = zeros(N,2,class_t);
err = zeros(N,4,class_t);
cnd = zeros(N,1,class_t);
d2 = zeros(N,4,class_t);
normsX = zeros(N,5,class_t);
rq2 = cell(d,1);
for j = 1:d
    rq2{j} = zeros(n,N,class_t);
end

p1 = zeros(n,1,class_t);
p2 = ones(n,1,class_t);
comp_eig = zeros(n,d,class_t);
restmp = zeros(N,d,class_t);
for j = 1:d
    norms(j) = norm(Aset{j});
end

realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(A{j}));
end

%% Experiments (N runs)
for k = 1:N
    Xs = cell(nb_LC, 1);
    Ys = cell(nb_LC, 1);
    lamus = cell(nb_LC, 1);
    mus = zeros(d, nb_LC); 
    % draw multiple linear combination of A
    for iter_lc = 1:nb_LC 
        % construct a normalized Gaussian vector mu
        if complex
            mu = randn(d,1,class_t) + 1i*randn(d,1,class_t);
        else
            mu = randn(d,1,class_t);
        end
        mu = mu/norm(mu);
    
        % construct a linear combination A(mu)
        mus(:,iter_lc) = mu;
        M = mu(1)*A{1};
        for j = 2:d
           M = M + mu(j)*A{j};
        end
    
        % add random Gaussian noise with Frobenius norm noise to M = A(mu)
        % size is relative to the norm of A(1)
        if noise>0
            if complex
                P = randn(n,class_t) + 1i*randn(n,class_t);
            else
                P = randn(n,class_t);
            end
            P = P/norm(P,'fro');
            M = M + noise/norm(A{1})*P;
        end
    
        [X,D,Y] = eig(M);  % compute eigenvectors of M = A(mu)
        
        lamus{iter_lc} = diag(D);
        s = 1./diag(Y'*X); % abs(s) are condition numbers of individual eigenvalues
        S = diag(s);       % factors for 2-sided RQ, as left and right eigenvectors are normalized  
        Y = Y*conj(S);      % ensures Y'*X = Id
        Xs{iter_lc} = X;
        Ys{iter_lc} = Y;
        
    end
    
    % combine eigenvectors found with matching function
    matching = greedyMatching(Xs);
    for t=2:nb_LC
        Xs{t} = Xs{t}(:, matching{t-1});
        Ys{t} = Ys{t}(:, matching{t-1});
    end

    % combine the matched eigenvectors
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
    
    cnd(k,1) = cond(X);
    s = 1./diag(Y'*X); % abs(s) are condition numbers of individual eigenvalues
    S = diag(s);       % factors for 2-sided RQ, as left and right eigenvectors are normalized  
    Y = Y*conj(S); 
    mu_avg = mean(mus,2);
    mu_avg = mu_avg / norm(mu_avg);
    M_avg = mu_avg(1)*A{1};
    for j = 2:d
       M_avg = M_avg + mu_avg(j)*A{j};
    end
    exact_lamu = exact * mu_avg;
    lamu = diag(Y' * M_avg * X);
    % find a reordering of the computed eigs that best matches the exact ones
    p = matchrows(exact_lamu,lamu,1,1,p1,p2);

    % we compute norms of the matrices D(mu), S(mu), T(mu) from the theory
    center_eig = sum(exact(clus,:),1)/clus_size; % 
    center_lamu = center_eig*mu;
    center_lamu_comp = sum(lamu(p(clus)))/clus_size; % 

    for ind = 1:clus_comp_size
        indeig = clus_comp(ind); 
        dif = norm(exact(indeig,:)-center_eig);
        stev_d(ind) = dif/abs(lamu(p(indeig))-center_lamu);
        stev_s(ind) = 1/abs(lamu(p(indeig))-center_lamu);
        stev_t(ind) = dif/abs(lamu(p(indeig))-center_lamu)^2;
        gap_s(ind) = 1/abs(lamu(p(indeig))-center_lamu_comp);
    end
    d2(k,:) = [norm(stev_d,Inf) norm(stev_s,Inf) norm(stev_t,Inf) norm(gap_s,Inf)];
    % we compute norms ||X1||, ||X2||, ||Y1||, ||Y2|| and ||D_k|| using diagonalization
    mX1A = X(:,p(clus));
    mX2A = X(:,p(clus_comp));
    mY1A = Y(:,p(clus));
    mY2A = Y(:,p(clus_comp));
    MA22 = mY2A'*M*mX2A - center_lamu_comp*eye(clus_comp_size);
    matnrm = [];
    for j = 1:d
        A22 = mY2A'*A{d}*mX2A - center_eig(d)*eye(clus_comp_size);
        matnrm(j) = norm(A22*inv(MA22));
    end

    normsX(k,:) = [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)];

    % 1-sided part -------------------------------------
    for j = 1:d
        comp_eig(:,j) = diag(X'*A{j}*X);
    end
    for j = 1:d
        rq1{j}(:,k) = comp_eig(p,j);
    end

    tmpres = compute_residual(Aset,comp_eig(p(clus),:),norms);
    tmpvec = tmpres(:,1).^2;
    for j = 2:d 
        tmpvec = tmpvec + tmpres(:,j).^2;
    end
    res(k,1) = norm(sqrt(tmpvec));
    err(k,1) = norm(comp_eig(p(clus),:)-exact(clus,:),'fro');
    err(k,2) = norm(comp_eig(p(clus(1)),:)-exact(clus(1),:),'fro');

    % 2-sided part -------------------------------------
    for j = 1:d
        comp_eig(:,j) = diag(Y'*A{j}*X);
    end
    
    for j = 1:d
        rq2{j}(:,k) = comp_eig(p,j);
    end

    warning on
    % residuals for all matrices of the family (n x d matrix)
    tmpres = compute_residual(Aset,comp_eig(p(clus),:),norms);
    tmpvec = tmpres(:,1).^2;
    for j = 2:d 
        tmpvec = tmpvec + tmpres(:,j).^2;
    end
    res(k,2) = norm(sqrt(tmpvec));
    err(k,3) = norm(comp_eig(p(clus),:)-exact(clus,:),'fro');
    err(k,4) = norm(comp_eig(p(clus(1)),:)-exact(clus(1),:),'fro');

end