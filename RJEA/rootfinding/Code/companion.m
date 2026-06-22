function [C] = companion(poly)
    % COMPANION constructs the companion matrix of a univariate polyomial.
    %
    %   [C] = COMPANION(p) constructs the companion matrix of a univariate
    %   polyomial in its matrix representation p.
    %
    %   [C] = COMPANION(P) constructs the companion matrix of a univariate
    %   polyomial in its cell representation P.
    %
    %   [C] = COMPANION(system) constructs the companion matrix of a
    %   univariate polynomial in its system representation system.
    %
    %   Input/output variables:
    %   
    %       C: [double(d,d)] companion matrix.
    %       p: [double(d+1,2)] coefficient matrix of the univariate
    %       polynomial.
    %       P: [cell(1,1)] cell representation of the univariate
    %       polynomial.
    %       system: [struct] system representation of the univariate
    %       polynomial.
    %
    %   See also ROOTS.
    
    % MacaulayLab (2022) - Christof Vermeersch.

    % TODO: think about use of multiple representations.
    
    % Transform the univariate polynomial into its system representation:
    if iscell(poly)
        poly = systemstruct(poly);
    elseif ismatrix(poly)
        poly = systemstruct({poly});
    end

    % Check if the polynomial is univariate:
    if poly.n ~= 1
        error('The polynomial is not univariate!')
    end
    
    % Construct the companion matrix:
    matrix = [poly.coef{1} poly.supp{1}];
    p = expandedpoly(matrix);
    C = zeros(poly.dmax,poly.dmax);
    C(2:poly.dmax,1:poly.dmax-1) = eye(poly.dmax-1);
    C(:,poly.dmax) = -p(1:end-1)'./p(end);
end