function p = matchrows(A,B,mu1,mu2,p1,p2)

% MATCHROWS(A,B)  returns perturbation p that minimizes norm(A-B(p,:))

% This is the assigmnent problem that could be solved by the Hungarian
% algorithm. Here we use a random approach using two random projections to
% a line. We repeat the method until we get the same ordering from both
%
% Optional arguments:
% mu1 and mu2 are initial random vectors (since random generator is slow,
% this way we can this use the same projection for all examples)
% p1 and p2 are initial vectors with zeros and ones

% Bor Plestenjak 2024

class_t = superiorfloat(A,B);
[n,d] = size(A);

if n<6
    % for n<6 we check all possible permutations, otherwise this is too slow
    P = perms(1:n);
    PB = B(P',:);
    best = Inf;
    for j = 1:size(P,1)
        tmp = norm(A-B(P(j,:)',:));
        if tmp<best
            best = tmp;
            p = P(j,:)';
        end
    end
else

if nargin<3
    mu1 = randn(d,1,class_t)+1i*randn(d,1,class_t);
end
if nargin<4
    mu2 = randn(d,1,class_t)+1i*randn(d,1,class_t);
end
if nargin<5
    p1 = zeros(n,1,class_t);
end
if nargin<6
    p2 = ones(n,1,class_t);
end

runs = 0;
while norm(p1-p2)>0 && runs<3
    runs = runs + 1;
    if runs>1
        % if we did not get the same order for both projections we repeat the
        % test with two new projections
        mu1 = randn(d,1,class_t)+1i*randn(d,1,class_t);
        mu2 = randn(d,1,class_t)+1i*randn(d,1,class_t);
    end
    pA1 = A*mu1;
    pB1 = B*mu1;
    pA2 = A*mu2;
    pB2 = B*mu2;
    for row = 1:n
        % we go through rows in a loop and assign each row to the nearest
        % nonassigned row from B
        [tmp,pos1] = min(abs(pB1-pA1(row)));
        [tmp,pos2] = min(abs(pB2-pA2(row)));
        p1(row) = pos1;
        pB1(pos1) = Inf;
        p2(row) = pos2;
        pB2(pos2) = Inf;
    end
end

if norm(p1-p2)>0 
    P = perms(1:n);
    PB = B(P',:);
    best = Inf;
    for j = 1:size(P,1)
        tmp = norm(A-B(P(j,:)',:));
        if tmp<best
            best = tmp;
            p = P(j,:)';
        end
    end
else
    p = p1;
end

end