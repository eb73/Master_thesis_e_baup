function [v] = stdmonomials(Z,tol)
    % STDMONOMIALS identifies the standard monomials.
    %
    %   [v] = STDMONOMIALS(Z) finds the standard monomials in the
    %   numerical basis matrix Z of the null space via rank checks.
    %
    %   [v] = STDMONOMIALS(___,tol) uses a user-specified tolerance for the
    %   required rank checks.
    %
    %   Input/output variables:
    %
    %       v: [int(mb,1)] list with standard monomials in Z.
    %       Z: [double(q,mb)] numerical basis of the null space.
    %       tol: [double - optional] tolerance of the rank tests (default 
    %       = 1e-10).
    %
    %   See also GAP, GAPSTDMONOMIALS.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Process the optional parameter:
    if ~exist('tol','var')
        tol = 1e-10; 
    end
    
    % Identify the standard monomials:
    ma = size(Z,2);
    vlist = zeros(size(Z,1),1);
    oldrank = 0;
    for k = 1:size(Z,1)
        newrank = rank(Z(1:k,:),tol);
        if newrank > oldrank
            vlist(k) = 1;
        end
        oldrank = newrank;
        if newrank == ma
            break
        end
    end
    v = find(vlist);
end