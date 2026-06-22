function res = compute_residual(Aset,lambda,norms)

% computes relative residual for a set of commuting matrices and eigenvalues. 
% We compute values min(svd(A-lambda*I))/norm(A)

% Optional arguments:
% norms is a vector with norms of matrices in Aset, if we provide norms the evaluation is faster

% Bor Plestenjak 2024

class_t = superiorfloat(Aset{:});

m = size(lambda,1);
k = length(Aset);
n = size(Aset{1},1);

if nargin<3
    res = zeros(m,k,class_t);
end
if nargin<4
    for j = 1:k
        norms(j) = norm(Aset{j});
    end
end
for j = 1:k
    A = Aset{j};
    for ind = 1:m
        res(ind,j) = min(svd(A-lambda(ind,j)*eye(n,class_t)))/norms(j);
    end
end