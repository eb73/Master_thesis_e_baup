function [pos] = posgrlex(V)
    % POSGRLEX implements the graded lexicographic monomial ordering.
    % 
    %
    %   [pos] = POSGRLEX(V) calculates the position of the monomial in
    %   the monomial basis using the graded lexicographic ordering.
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
        d = sum(v);
        order = nbmonomials(d,nc);
        for l = 1:nc
            order = order - nbmonomials(d-sum(v(1:l))-1,nc-l);
        end
        pos(k) = order;
    end
end