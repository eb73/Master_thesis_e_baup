function [pos] = position(V,order)
    % POSITION determines the position of a monomial.
    %
    %   [pos] = POSITION(V) calculates the position of the monomial in
    %   the monomial basis using the graded inverse lexicographic ordering.
    %
    %   [___] = POSITION(___,order) uses a user-defined monomial ordering. 
    %
    %   Input arguments:
    %       V [int(k,n)]: monomial vector representation of the monomials
    %       (without the coefficients).
    %       order [function - optional]: user-defined monomial ordering
    %       (default = @posgrinvlex).
    %
    %   Output argument:
    %       pos [int(k,1)]: position in the monomial basis.   
    %
    %   See also MONOMIALSMATRIX.
            
    % MacaulayLab (2023) - Christof Vermeersch.
    
    % Set the default monomial ordering:
    if ~exist('order','var') 
        order = @posgrinvlex;
    end
    
    pos = order(V);
end