function [shifted,coefficients] = basischeb(original,shift)
    % BASISCHEB is the multiplication property of the Chebyshev polynomial
    % basis.
    %
    %   [shifted,coefficients] = BASISCHEB(original,shift) multiplies the
    %   original monomial with the shift polynomial in the Chebyshev
    %   polynomial basis
    %
    %   Input arguments:
    %       original [doubles(1,n)]: support of the original monomial.
    %       shift [doubles(k,n)]: support of the shift polynomial.
    %
    %   Output arguments:
    %       shifted [doubles(k*l,n)]: support of the shifted polynomial.
    %       coefficients [double(l,1)]: coefficients of the multiplication
    %       property.
    %
    %   See also EVALCHEB.

    % MacaulayLab (2023) - Christof Vermeersch.
    
    shifted = original;
    for k = 1:size(shift,2)
        % Perform the shift in every variable sequentially. 
        shiftx = zeros(size(shift));
        shiftx(:,k) = shift(:,k);
        shifted = [shifted + shiftx; abs(shifted - shiftx)];
    end
    coefficients = (1/(2^size(shift,2)))*ones(2^size(shift,2),1);
end