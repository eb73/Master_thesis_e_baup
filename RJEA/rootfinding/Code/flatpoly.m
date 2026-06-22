function q = flatpoly(p)
    % FLATPOLY flattens a polynomial in its coefficient matrix
    % representation.
    %
    %   [q] = FLATPOLY(p) sums corresponding entries (rows) of the
    %   coefficient matrix representation.
    %
    %   Input/output variables:
    %
    %       q: [double(l,n+1)] flatened coefficient matrix.
    %       p: [double(k,n+1)] coefficient matrix.

    % MacaulayLab (2022) - Christof Vermeersch.
    
    [~, unq, corr] = unique(p(:,2:end),'rows');
    idxduplicates = setdiff(1:size(p,1), unq);
    q = p(unq,:);
    for k = 1:length(idxduplicates)
        q(corr(idxduplicates(k)),1) = q(corr(idxduplicates(k)),1) ...
             + p(idxduplicates(k),1);
    end
end