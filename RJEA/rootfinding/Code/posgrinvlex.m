function [pos] = posgrinvlex(V)
    % POSGRINVLEX implements the graded inverse lexicographic monomial
    % ordering.
    %
    %   [pos] = POSGRINVLEX(V) calculates the position of the monomial in
    %   the monomial basis using the graded inverse lexicographic ordering.
    %
    %   Input argument:
    %       V [int(k,n)]: monomial vector representation of the monomials
    %       (without the coefficients).
    %
    %   Output argument:
    %       pos [int(k,1)]: position in the monomial basis.   
    %
    %   See also POSITION.
            
    % MacaulayLab (2023) - Christof Vermeersch.
    
    % TODO: vectorize code.
    
    [nr,nc] = size(V);
    pos = zeros(1,nr);
    for k = 1:nr
        v = V(k,:);
        order = 0;
        for l = 1:nc-1
            order = order + nbmonomials(sum(v(l:nc),2)-1,nc-l+1);
        end
        pos(k) = order + v(nc) + 1;
    end
end