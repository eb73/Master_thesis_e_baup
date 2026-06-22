function [q] = contractedpoly(p,d,n,tol)
    % CONTRACTEDPOLY contracts a polynomial into its corresponding
    % coefficient matrix.
    %
    %   [q] = CONTRACTEDPOLY(p,d,n) contracts a polynomial p (with maximum
    %   total degree d and n variables) given by its vector representation 
    %   (in the degree negative lexicographic ordering) into its 
    %   corresponding coefficient matrix q.
    %
    %   [q] = CONTRACTEDPOLY(___,tol) only considers coefficients larger
    %   than the user-specified tolerance.
    %
    %   Input/output variables:
    %
    %       q: [double(k,n+1)] corresponding coefficient matrix.
    %       p: [double(1,l)] row vector.
    %       d: [int] degree of the polynomial.
    %       n: [int] number of variables.
    %       tol: [double - optional] tolerance for a non-zero coefficient
    %       (default = 1e-10).
    %
    %   See also EXPANDEDPOLY.
    
    % MacaulayLab (2022) - Christof Vermeersch.  

    % Process the optional parameter:
    if ~exist('tol','var') 
        tol = 1e-10;
    end
    
    % Contract the polynomial:
    mask = find(abs(p) > tol);
    K = monomialsmatrix(d,n);
    q = [p(mask)' K(mask,:)];    
end