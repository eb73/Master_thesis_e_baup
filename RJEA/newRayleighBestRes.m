function [rq1, rq2, res, err, cnd, normsX, d2] = newRayleighBestRes(A,N,exact,clus,complex,noise,Aset,nb_LC)
% Computes the approximate joint eigenvalues with RJEA+best residual(k=nb_LC)
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
%  - rq2: same as above when 2-sided RQ is used
%  - res: vector of size N x 2 with residuals
%         v(k) = norm((A{k}-rho(j,k)*I)*X(:,j))/norm(A{k}), 
%         where rho(j,k) is an eigenvalue approximation (2RQ and 1RQ)
%  - err: N x 4 matrix with 2-norm of the difference between the computed 
%         joint eigenvalue and the exact eigenvalue (2RQ and 1RQ) (first
%         for cluster and then for the first index from the cluster)
%  - cnd: vector of size N with condition numbers of eigenvector matrix X
%  - d2: values d(mu), s(mu), t(mu) and gap
%  - normsX : values [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)] from diagonalization

class_t = superiorfloat(A{:});
n = size(A{1},1);     % size of commuting matrices 
d = length(A);        % number of commuting matrices in set A 

res = zeros(N,2,class_t);
err = zeros(N,4,class_t);
cnd = zeros(N,1,class_t);
d2 = zeros(N,4,class_t);
normsX = zeros(N,5,class_t);
rq1 = cell(d,1);
rq2 = cell(d,1);
for j = 1:d
    rq1{j} = zeros(n,N,class_t);
    rq2{j} = zeros(n,N,class_t);
end
for k = 1:N
    best_res1 = Inf;
    best_res2 = Inf;
    for iter_lc = 1:nb_LC
        [rq2_tmp, rq1_tmp, res_tmp, err_tmp, cnd_tmp, d2_tmp, normsX_tmp, ~] = testRayleigh(A,1,exact,clus,complex,noise,0,Aset,0);
        res_tmp1 = res_tmp(:, 2);
        res_tmp2 = res_tmp(:, 1);
        % best residual for RQ1
        if res_tmp1 < best_res1
            best_res1 = res_tmp1;
            for j = 1:d
                rq1{j}(:,k) = rq1_tmp{j}(:,1);
            end
            res(k,1) = res_tmp1;
            err(k,[1,2]) = err_tmp(:, [2,4]);
            cnd(k,:) = cnd_tmp;
            d2(k,:) = d2_tmp;
            normsX(k,:) = normsX_tmp;
        end
        % best residual for RQ2
        if res_tmp2 < best_res2
            best_res2 = res_tmp2;
            for j = 1:d
                rq2{j}(:,k) = rq2_tmp{j}(:,1);
            end
            res(k,2) = res_tmp2;
            err(k,[3,4]) = err_tmp(:, [1,3]);
        end

    end
end