function [pos] = posgrevlex(V)
    % POSGREVLEX implements the graded reverse lexicographic monomial 
    % ordering.
    % 
    %
    %   [pos] = POSGREVLEX(V) calculates the position of the monomial in
    %   the monomial basis using the graded reverse lexicographic ordering.
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
        v = fliplr(V(k,:));
        d = sum(v);
        order = nbmonomials(d,nc);
        for l = 1:nc
            order = order - nbmonomials(d-sum(v(1:l))-1,nc-l);
        end
        pos(k) = order;
    end
end