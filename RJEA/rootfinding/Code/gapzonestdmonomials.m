function [dgap,ma,gapsize] = gapzonestdmonomials(c,dend,n,blocksize)
    % GAPZONESTDMONOMIALS identifies the gap via the standard monomials.
    %
    %   [dgap,ma,gapsize] = GAPZONESTDMONOMIALS(c,dend,n) returns the
    %   degree of the first degree block of the gap zone, the number of
    %   affine solutions, and the size of the gap zone based on the list of 
    %   identified standard monomials c.
    %
    %   [___] = GAPZONESTDMONOMIALS(___,blocksize) considers the number of 
    %   rows in the block rows (cf. the null space of block Macaulay 
    %   matrix).
    %
    %   Input/output variables:
    %  
    %       dgap: [int] degree of first degree block of the gap zone.
    %       ma: [int] number of affine solutions.
    %       gapsize: [int] number of degree blocks in the gap zone.
    %       c: [int(k,1)] list with positions of the identified standard 
    %       monomials (in the degree negative lexicographic ordering).  
    %       dend: [int] highest total degree of the identified standard 
    %       monomials.
    %       n: [int] number of variables.
    %       blocksize: [int - optional] number of rows in the block rows 
    %       (necessary for the block Macaulay matrix - default = 1).
    %
    %   See also GAP, GAPZONE, STDMONOMIALS.
    %
    %   Note: GAPZONESTDMONOMIALS does not throw an error when there is no 
    %   gap zone.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Process the optional parameter:
    if ~exist('blocksize','var') 
        blocksize = 1;
    end
    
    % Determine gap degree,number of affine solutions, and size of the gap:
    ma = 0;
    gapsize = 0;
    dgap = NaN;
    for k = 0:dend
        nm = blocksize*nbmonomials(k,n);
        r = sum(c<=nm);
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