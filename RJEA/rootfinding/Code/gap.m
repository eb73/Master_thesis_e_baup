function [dgap,ma] = gap(Z,dend,n,blocksize,tol)
    % GAP finds the gap zone in the null space via a basis matrix.
    %
    %   [dgap,ma] = GAP(Z,dend,n) returns the first degree dgap of the gap
    %   zone in a basis matrix of the null space Z (and the number of affine
    %   solutions ma).
    %
    %   [___] = GAP(___,blocksize) considers the number of rows in the
    %   block rows (cf. the null space of block Macaulay matrix).
    %
    %   [___] = GAP(___,tol) uses a user-specified tolerance for the
    %   required rank checks.
    %
    %   Input/output variables:
    %  
    %       dgap: [int] first degree of the gap zone.
    %       ma: [ma] number of affine solutions.
    %       Z: [double(q,mb)] basis matrix of the null space of a (block) 
    %       Macaulay matrix.
    %       dend: [int] highest total degree of basis matrix.
    %       n: [int] number of variables.
    %       blocksize: [int - optional] number of rows in the block rows 
    %       (necessary for the block Macaulay matrix - default = 1).
    %       tol: [double - optional] tolerance of the rank tests (default 
    %       = 1e-10).
    %
    %   See also GAPSTDMONOMIALS, STDMONOMIALS.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Process the optional parameters:
    if ~exist('tol','var') 
        tol = 1e-10;
    end
    if ~exist('blocksize','var') 
        blocksize = 1;
    end
    
    % Determine gap degree and number of affine solutions:
    ma = 0;
    for k = 0:dend
        nm = blocksize*nbmonomials(k,n);
        r = rank(Z(1:nm,:),tol);
        if ma == r
            dgap = k;
            break
        else
            ma = r;
        end
    end
    
    % Throw an error if there is no gap:
    if ~exist('dgap','var') 
        error('No gap found! Maybe increase the degree?')
    end
end