function [rq1, rq2, res1, err1, res2, err2, cnd, d2, normsX] = newRayleighOpt(A,Aset,N,complex,param)

% NEWRAYLEIGHOPT(A,Aset,N,complex,param) computes the joint eigenvalues approximation of A for N runs, 
% using the method and setting specified in param.
% 
% Input:
%  - A: a set of d commuting n x n matrices (usually perturbed matrices from Aset)
%  - Aset: original matrices used for computing residuals
%  - N: number of runs
%  - complex: (1) use complex random coefficients, 0: use real
%  - param: struct containing setting (method, nb_trials, rq_type, etc)
%  
% Output:
%  - rq1: cell of size d with n X N matrices, where (i,j)-th element of k-th
%        cell is 1-sided RQ approximation for exact(i,k) in run j.
%        if rq_type='two-sided', rq1 is empty
%  - rq2: same as above but when 2-sided RQ is used.
%        if rq_type='one-sided', rq2 is empty
%  - res1: vector of size N x 1 with residuals for RQ1 
%          if rq_type='two-sided', contains all zeros (and is unused in main)
%  - res2: vector of size N x 1 with residuals for RQ2
%          if rq_type='one-sided', contains all zeros 
%  - err1: N x 2 matrix with 2-norm of the difference between the computed 
%         joint eigenvalue from RQ1 and the exact eigenvalue (first
%         for cluster and then for the first index from the cluster)
%         if rq_type='two-sided', contains all zeros
%  - err2: same as above with RQ2 errors
%         if rq_type='one-sided', contains all zeros
%  - cnd: vector of size N with condition numbers of eigenvector matrix X
%  - d2: values d(mu), s(mu), t(mu) and gap
%  - normsX: values [norm(mX1A) norm(mX2A) norm(mY1A) norm(mY2A) norm(matnrm)] from diagonalization

%% Setup

class_t = superiorfloat(A{:});
d = length(A); % number of commuting matrices

exact = param.exact_eig;
clus = param.index;
noise=  param.noise_for_Amu;
nb_LC = param.nb_trials; 
rq_type = param.rq_type;
method = param.method;

realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(A{j}));
end

rq1 = {};
rq2 = {};
res1 = zeros(N,1,class_t);
err1 = zeros(N,2,class_t);
res2 = zeros(N,1,class_t);
err2 = zeros(N,2,class_t);

%% Redirect experiment to the method in param

switch method
    case 'best_res'
        [rq1, rq2, res, err, cnd, normsX, d2] = newRayleighBestRes(A,N,exact,clus,complex,noise,Aset,nb_LC);
        res1 = res(:, 1);
        res2 = res(:, 2);
        err1 = err(:, [1,2]);
        err2 = err(:, [3,4]);
    case 'mean'
        [rq1, rq2, res, err, cnd, d2, normsX] = newRayleighMean(A,N,exact,clus,complex,noise,Aset,nb_LC);
        res1 = res(:, 1);
        res2 = res(:, 2);
        err1 = err(:, [1,2]);
        err2 = err(:, [3,4]);
    case 'smaller_subspace'
        % for this method, need to specify type of rayleigh quotient
        if ismember(rq_type,  {'one-sided', 'both'})
            [rq1, res1, err1, cnd1, d2, normsX] = newRayleighSmallerSubspace_rq1(A,N,exact,clus,complex,noise,Aset, nb_LC);
            cnd = cnd1;
        end
        if ismember(rq_type, {'two-sided', 'both'})
            [rq2, res2, err2, cnd2, d2, normsX] = newRayleighSmallerSubspace_rq2(A,N,exact,clus,complex,noise,Aset, nb_LC);
            cnd = cnd2;
        end
        % the condition number in this case aren't the same for rq1/rq2 so keep track of both
        if strcmp(rq_type, 'both')
            cnd = [cnd1, cnd2];
        end
    case 'ref'
        [~, rq1, res, err, cnd, d2, normsX, ~] = testRayleigh(A,N,exact,clus,complex,noise,0,Aset,0);
        % separate RQ1 and RQ2 res/err (the indices in the original code testRayleigh are not the same)
        res1 = res(:, 2);
        res2 = res(:, 1);
        err1 = err(:, [2,4]);
        err2 = err(:, [1,3]);
    otherwise
        error(['Method ', method, ' is not compatible, please choose one of {"best_res", "mean", "smaller_subspace"}'])
end

