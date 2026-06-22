function [matchings] = greedyMatching(Xs)
% Performs a greedy matching of the columns of X1 = Xs{1} with the columns of the
% other matrices Xs{t}. This is used to match the eigenvectors in different eigenbases
% with the same eigenvalue.
%
% Input: 
%  - Xs (cell): contains T n x n matrices (that represent different
%  eigenbases, typically obtained through different trials in RJEA)

% Output:
%  - matchings: cell of T-1 vectors of length n. matchings{t} corresponds to 
%  matching the columns of X1 with the cols of Xs{t+1} (i.e. to match eigenvectors 
%  in different eigenbases associated with the same eigenvalue).

T = numel(Xs);

if T < 2
    error('Cell has size smaller than 2 (to do a matching, you must have at least 2 eigenbases).')
end
size1 = size(Xs{1});
for t=2:T
    if ~isequal(size(Xs{t}), size1)
        error('Cells must have same sizes, there is a  mismatch between cell 1 and ' + str(t))
    end
end


n = length(Xs{1});
matchings = cell(T-1, 1);
X1t = Xs{1}';

for t = 2:T
    matchings{t-1} = zeros(n,1);
    S = abs(X1t*Xs{t});
    for iter =1:n
        max_val = max(max(S));
        [max_row,max_col] = find(S==max_val);
        i = max_row(1);
        j = max_col(1);
        matchings{t-1}(i) = j; % match col i of X1 with col j of Xt
        S(i,:) = -Inf;
        S(:,j) = -Inf;
    end
end

end