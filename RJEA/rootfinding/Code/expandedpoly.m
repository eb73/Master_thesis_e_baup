function [q] = expandedpoly(p)
    % EXPANDEDPOLY expands a polynomial in its row vector representation.
    %
    %   [q] = EXPANDEDPOLY(p) expands a polynomial p given by its 
    %   coefficient matrix p into the corresponding row vector q (in the 
    %   degree negative lexicographic ordering).
    %
    %   Input/output variables:
    %
    %       q: [double(1,l)] corresponding row vector.
    %       p: [double(k,n+1)] coefficient matrix.
    %
    %   See also CONTRACTEDPOLY.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Create an empty expanded polynomial of appropriate length:
    n = size(p,2) - 1;
    d = max(sum(p(:,2:end),2));    
    q = zeros(1,nbmonomials(d,n));
    
    % Add non-zero coefficients to the expanded polynomial:
    coef = p(:,1);
    pos = position(p(:,2:end));
    q(pos) = coef;
end