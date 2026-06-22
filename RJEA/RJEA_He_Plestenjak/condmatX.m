function A = condmatX(n,kappa,class_t)

% A = CONDMATX(n,kappa) random n x n real matrix with condition number kappa
% and normalized vectors (that could be used as an eigenvector matrix with
% a prescribed condition number)

% Bor Plestenjak 2024

if nargin<3
    class_t = 'double';
end

if n>1000
    error('n is too large')
end

if kappa<1
    error('Condition number can not be below 1')
end

[Q1,R] = qr(randn(n,class_t));
[Q2,R] = qr(randn(n,class_t));

function y = f(z)
   D = diag(z.^(-(0:n-1)/(n-1)));
   X = Q1*D*Q2;
   X1 = X/diag(diag(sqrt(X'*X)));
   y = cond(X1) - kappa;
end

kappacor = fzero(@f,kappa);

D = diag(kappacor.^(-(0:n-1)/(n-1)));
X = Q1*D*Q2;
A = X/diag(diag(sqrt(X'*X)));
end