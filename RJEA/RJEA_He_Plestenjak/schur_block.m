function [X1,X2,Y1,Y2] = schur_block(M,exact,indices)

% Returns X = [X1 X2] such that inv(X)*M*X is block diagonal with eigenvalues
% exact with given indices in block (1,1) and remaining eigenvalues in
% block (2,2)

% We use ordered Schur form and transformation P = [I X; 0 I], where X is a
% solution of a Sylvester equation

% Bor Plestenjak 2024

m = length(indices);
n = size(M,1);
[Q,R] = schur(M,'complex');
r = diag(R);
p = matchrows(exact(:),r(:));
ordsch = zeros(1,n);
ordsch(p(indices)) = 1;
[Q1,R1] = ordschur(Q,R,ordsch);
X = lyap(R1(1:m,1:m),-R1(m+1:n,m+1:n),R1(1:m,m+1:n));
P = eye(n);
P(1:m,m+1:n) = X;
Z = Q1*P;
warning off
B = inv(P)*R1*P;
warning on
X1 = Z(:,1:m);
X2 = Z(:,m+1:n);
for k = 1:size(X2,2)
    X2(:,k) = X2(:,k)/norm(X2(:,k));
end
Y = inv([X1 X2]');
Y1 = Y(:,1:m);
Y2 = Y(:,m+1:n);




