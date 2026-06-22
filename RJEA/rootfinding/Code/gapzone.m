function [dgap,ma,gapsize] = gapzone(Z,dend,n,tol,recr,blocksize)
    % GAPZONE identifies the gap.
    %
    %   [dgap,ma,gapsize] = GAPZONE(Z,dend,n) returns the degree of the
    %   first degree block of the gap zone, the number of affine solutions,
    %   and the size of the gap zone.
    %
    %   [___] = GAPZONE(___,tol) uses a user-specified tolerance for
    %   the required rank checks.
    %
    %   [___] = GAPZONE(___,tol,blocksize) considers the number of rows 
    %   in the block rows (cf. the null space of block Macaulay matrix).
    %
    %   Input/output variables:
    %  
    %       dgap: [int] degree of first degree block of the gap zone.
    %       ma: [int] number of affine solutions.
    %       gapsize: [int] number of degree blocks in the gap zone.
    %       Z: [double(q,mb)] basis matrix of the null space of a (block) 
    %       Macaulay matrix.
    %       dend: [int] highest total degree of basis matrix.
    %       n: [int] number of variables.
    %       tol: [double - optional] tolerance of the rank tests (default 
    %       = 1e-10).
    %       blocksize: [int - optional] number of rows in the block rows 
    %       (necessary for the block Macaulay matrix - default = 1).
    %
    %   See also GAP, GAPSTDMONOMIALS, STDMONOMIALS, GAPZONESTDMONOMIALS.
    %
    %   Note: GAPZONE does not throw an error when there is no gap zone.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Process the optional parameters:
    if ~exist('tol','var') 
        tol = 1e-10;
    end
    if ~exist('recr','var') 
        recr = false;
    end
    if ~exist('blocksize','var') 
        blocksize = 1;
    end
    
    % Determine gap degree,number of affine solutions, and size of the gap:
    gapsize = 0;
    dgap = NaN;
    if recr
        q = size(Z,2);
        nm = blocksize*nbmonomials(0,n);
        Y = null(Z(1:nm,:),tol);
        ma = q - size(Y,2); 
        for k = 1:dend
            nm1 = blocksize*nbmonomials(k-1,n)+1;
            nm2 = blocksize*nbmonomials(k,n);
            Y = nullrecrrow(Y,Z(nm1:nm2,:));
            r = q - size(Y,2);
            if ma == r
                gapsize = gapsize + 1;
                dgap = k - gapsize + 1;
            else
                if ~isnan(dgap)
                    break
                else
                    ma = r;
                end
            end
        end
    else
        ma = 0;
        for k = 0:dend
            nm = blocksize*nbmonomials(k,n);
            r = rank(Z(1:nm,:),tol);
            if ma == r
                gapsize = gapsize + 1;
                dgap = k - gapsize + 1;
            else
                if ~isnan(dgap)
                    break
                else
                    ma = r;
                end
            end
        end
    end
end