function [s] = nbmonomials(d,n)
    % NBMONOMIALS calculates the number of monomials.
    %
    %   [s] = NBMONOMIALS(d,n) calculates the number of monomials s in the 
    %   monomial basis for n variables and maximum total degree d.
    %
    %   Input/output variables:
    %
    %       s: [int] number of monomials in the monomial basis.
    %       d: [int] maximum total degree of the monomials.
    %       n: [int] number of variables of the monomials.
    %
    %   See also MONOMIALSMATRIX.
    
    % MacaulayLab (2022) - Christof Vermeersch.
    
    if d < 0
        s = 0; % Avoid errors with negative degrees.
    else
        s = nchoosek(d+n,n); % Number of monomials in the monomial basis.
    end
end