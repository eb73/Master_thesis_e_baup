function [shifted,coefficients] = basismon(original,shift)
    % BASISMON is the multiplication property of the standard monomial
    % basis.
    %
    %   [shifted,coefficients] = BASISMON(original,shift) multiplies the
    %   original monomial with the shift polynomial in the standard 
    %   monomial basis.
    %
    %   Input arguments:
    %       original [doubles(1,n)]: support of the original monomial.
    %       shift [doubles(k,n)]: support of the shift polynomial.
    %
    %   Output arguments:
    %       shifted [doubles(k,n)]: support of the shifted polynomial.
    %       coefficients [double]: coefficients of the multiplication
    %       property.
    %
    %   See also EVALMON.

    % MacaulayLab (2023) - Christof Vermeersch.

    shifted = original + shift;
    coefficients = 1;
end