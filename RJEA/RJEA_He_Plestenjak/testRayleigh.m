function [rq2, rq1, res, err, cnd, d2, normsX, SA] = testRayleigh(A,N,exact,clus,complex,noise,no_be,Aset,epslevel)

% [res, err, rq, cnd, cndev, res1, err1, rq1] = testRayleigh(A,N,exact,clus,complex,noise,no_be,Aset,epslevel)
% performs N runs of the algorithm to compute eigenvalues of a family of
% commuting matrices A using eigenvectors of a random linear combination A(mu)
% and returns norms of the residuals, errors, condition number, ...
%
% Input:
%  - A: a set of d commuting n x n matrices (usually perturbet matrices from Aset)
%  - N: number of runs
%  - exact: n x d matrix with exact joint eigenvalues of A for comparison
%  - clus: index of eigenvalue or vector of indices of a cluster
%  - complex: (1) use complex random coefficients, 0: use real
%  - noise: (0) amount of Gaussian noise added to A(mu) before computing eigenvector matrix X
%  - no_be: (0) 
%       0: use double preciosion (the precision of input matrices)
%       1: we use multiple precision to compute Rayleigh quotients without round-off error
%       2: we compute eigenvectors of A(mu) and Rayleigh quotients in multiple precision 
%       options 1 and 2 require Advanpix Multiprecision Computing Toolbox
%       (you can use slower VPA by setting the flag use_advanpix below to 0)
%  - Aset: (A) original matrices used for computing residuals
%  - epslevel: (0) noise added to matrices A (different perturbation in each run)
%  
% Output:
%  - rq2: cell of size d with n X N matrices, where (i,j)-th element of k-th
%        cell is 2-sided RQ approximation for exact(i,k) in run j
%  - rq1: same as rq when 1-sided RQ is used
%  - res: vector of size N x 2 with residuals
%         v(k) = norm((A{k}-rho(j,k)*I)*X(:,j))/norm(A{k}), 
%         where rho(j,k) is an eigenvalue approximation (2RQ and 1RQ)
%  - err: N x 4 matrix with 2-norm of the difference between the computed 
%         joint eigenvalue and the exact eigenvalue (2RQ and 1RQ) (first
%         for cluster and then for the first index from the cluster)
%  - cnd: vector of size N with condition numbers of eigenvector matrix X
%  - d2: values d(mu), s(mu), t(mu) and gap
%  - normsX : values [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)] from diagonalization
%  - SA _ values [norm(XX1) norm(XX2) norm(YY1) norm(YY2) norm(matnrm)] from Schur norm

% Bor Plestenjak 2024

class_t = superiorfloat(A{:});

% Change this to 0 to use slower vpa instead of Multiprecision Computing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_advanpix = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default values 
if nargin<5, complex = 1; end 
if nargin<6, noise = 0;   end 
if nargin<7, no_be = 0; end
if nargin<8, Aset = A; end
if nargin<9, epslevel=0; end

n = size(A{1},1);     % size of commuting matrices 
d = length(A);        % number of commuting matrices in set A 

% we compute complements if indices in clusters A and B
clus_comp = 1:n;
clus_comp(clus) = 0;
clus_comp = find(clus_comp~=0);
clus_comp_size = length(clus_comp);
clus_size = length(clus);

res = zeros(N,2,class_t);
err = zeros(N,4,class_t);
cnd = zeros(N,1,class_t);
res1 = zeros(N,n,class_t);
err1 = zeros(N,n,class_t);
d2 = zeros(N,4,class_t);
normsX = zeros(N,5,class_t);
SA = zeros(N,5,class_t);
rq2 = cell(d,1);
rq1 = cell(d,1);
for j = 1:d
    rq1{j} = zeros(n,N,class_t);
    rq2{j} = zeros(n,N,class_t);
end

testmu1 = randn(d,1,class_t) + 1i*randn(d,1,class_t); % random vector used for matching computed to exact eigenvalues
testmu2 = randn(d,1,class_t) + 1i*randn(d,1,class_t); % random vector used for matching computed to exact eigenvalues

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

% Experiments (N runs)
for k = 1:N
    if epslevel>0
        PA = cell(d,1);
        for j = 1:d 
            if realfamily
                PA{j} = randn(n);
            else
                PA{j} = randn(n) + 1i*randn(n);
            end
            PA{j} = PA{j}/norm(PA{j},'fro');
        end
        A = Aset;
        for j = 1:d 
            A{j} = A{j} + epslevel/sqrt(d)*PA{j};
        end
    end
    % construct a normalized Gaussian vector mu
    if complex
        mu = randn(d,1,class_t) + 1i*randn(d,1,class_t);
    else
        mu = randn(d,1,class_t);
    end
    mu = mu/norm(mu);

    % construct a linear combination A(mu)
    exact_lamu = mu(1)*exact(:,1);
    M = mu(1)*A{1};
    for j = 2:d
       M = M + mu(j)*A{j};
       exact_lamu = exact_lamu + mu(j)*exact(:,j);
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

    if no_be==2
        if use_advanpix
            [X,D,Y] = eig(mp(M));  % compute eigenvectors of M = A(mu)
        else
            [X,D] = eig(vpa(M));  % compute eigenvectors of M = A(mu)
            Y = inv(X');
            for col = 1:n
                Y(:,col) = Y(:,col)/norm(Y(:,col));
            end
        end
        lamu = double(diag(D)); % eigenvalue of A(mu)
    else
        [X,D,Y] = eig(M);  % compute eigenvectors of M = A(mu)
        lamu = diag(D); % eigenvalue of A(mu)
    end
    % abs(s) are condition numbers of individual eigenvalues
    s = 1./diag(Y'*X); % abs(s) are condition numbers of individual eigenvalues
    S = diag(s);       % factors for 2-sided RQ, as left and right eigenvectors are normalized   
    Y1 = Y*conj(S);    % left eigenvectors scaled so that Y1'*X = Id
    cnd(k,1) = cond(X);

    % find a reordering of the computed eigs that best matches the exact ones
    p = matchrows(exact_lamu,lamu,1,1,p1,p2);
    % 2-sided part -------------------------------------
    % compute eigenvales from RQ
    for j = 1:d
        if no_be==2
            if use_advanpix
                comp_eig(:,j) = double(diag(S*Y'*mp(A{j})*X));
            else
                comp_eig(:,j) = double(diag(S*Y'*vpa(A{j})*X));
            end
        elseif no_be==1
            if use_advanpix
                comp_eig(:,j) = double(diag(mp(S)*mp(Y)'*mp(A{j})*mp(X)));
            else
                comp_eig(:,j) = double(diag(vpa(S)*vpa(Y)'*vpa(A{j})*vpa(X)));
            end
        else
             comp_eig(:,j) = diag(S*Y'*A{j}*X);
        end
    end
    
    for j = 1:d
        rq2{j}(:,k) = comp_eig(p,j);
    end

    % we compute norms of the matrices D(mu), S(mu), T(mu) from the theory
    center_eig = sum(exact(clus,:),1)/clus_size; % 
    center_eig_kron = kron(center_eig,ones(clus_size,1));
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
    mY1A = Y1(:,p(clus));
    mY2A = Y1(:,p(clus_comp));
    Y1A(k) = norm(mY1A);
    Y2A(k) = norm(mY2A);
    MA22 = mY2A'*M*mX2A - center_lamu_comp*eye(clus_comp_size);
    matnrm = [];
    for j = 1:d
        A22 = mY2A'*A{d}*mX2A - center_eig(d)*eye(clus_comp_size);
        matnrm(j) = norm(A22*inv(MA22));
    end

    normsX(k,:) = [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)];
    
    % we compute norms ||X1||, ||X2||, ||Y1||, ||Y2|| and norm(D) from Schur form
    [XX1,XX2,YY1,YY2] = schur_block(M,exact_lamu,clus);
    MA22 = YY2'*M*XX2;
    MA11 = YY1'*M*XX1;
    matnrm = [];
    warning off
    for j = 1:d
        A22 = YY2'*A{j}*XX2 - center_eig(j)*eye(clus_comp_size);
        matnrm(j) = norm(kron(eye(n-clus_comp_size),A22)*inv(kron(eye(n-clus_comp_size),MA22)-kron(MA11.',eye(clus_comp_size))));
    end
    warning on
    SA(k,:) = [norm(XX1) norm(XX2) norm(YY1) norm(YY2) norm(matnrm)];

    % residuals for all matrices of the family (n x d matrix)
    tmpres = compute_residual(Aset,comp_eig(p(clus),:),norms);
    tmpvec = tmpres(:,1).^2;
    for j = 2:d 
        tmpvec = tmpvec + tmpres(:,j).^2;
    end
    res(k,1) = norm(sqrt(tmpvec));
    err(k,3) = norm(comp_eig(p(clus(1)),:)-exact(clus(1),:),'fro');
    err(k,1) = norm(comp_eig(p(clus),:)-exact(clus,:),'fro');
      
    % 1-sided part -------------------------------------
    for j = 1:d
        if no_be==2
            if use_advanpix
                comp_eig(:,j) = double(diag(X'*mp(A{j})*X));
            else
                comp_eig(:,j) = double(diag(X'*vpa(A{j})*X));
            end
        elseif no_be==1
            if use_advanpix
                comp_eig(:,j) = double(diag(mp(X)'*mp(A{j})*mp(X)));
            else
                comp_eig(:,j) = double(diag(vpa(X)'*vpa(A{j})*vpa(X)));
            end
        else
            comp_eig(:,j) = diag(X'*A{j}*X);
        end
    end
    for j = 1:d
        rq1{j}(:,k) = comp_eig(p,j);
    end

    tmpres = compute_residual(Aset,comp_eig(p(clus),:),norms);
    tmpvec = tmpres(:,1).^2;
    for j = 2:d 
        tmpvec = tmpvec + tmpres(:,j).^2;
    end
    res(k,2) = norm(sqrt(tmpvec));
    err(k,4) = norm(comp_eig(p(clus(1)),:)-exact(clus(1),:),'fro');
    err(k,2) = norm(comp_eig(p(clus),:)-exact(clus,:),'fro');

end