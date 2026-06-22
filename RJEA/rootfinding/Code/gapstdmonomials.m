function [dgap,ma] = gapstdmonomials(c,dend,n,blocksize)
    % GAPSTDMONOMIALS finds the gap zone in the null space via the standard 
    % monomials.
    %
    %   [dgap,ma] = GAPSTDMONOMIALS(c,dend,n) returns the first degree dgap
    %   of the gap zone in a basis matrix of the null space (and the number 
    %   of affine solutions ma) based on the list of identified standard 
    %   monomials c.
    %
    %   [dgap,ma] = GAPSTDMONOMIALS(___,blocksize) considers the number of 
    %   rows in the block rows (cf. the null space of block Macaulay 
    %   matrix).
    %
    %   Input/output variables:
    %  
    %       dgap: [int] first degree of the gap zone.
    %       ma: [ma] number of affine solutions.
    %       c: [int(k,1)] list with positions of the identified standard 
    %       monomials (in the degree negative lexicographic ordering).  
    %       dend: [int] highest total degree of the identified standard 
    %       monomials.
    %       n: [int] number of variables.
    %       blocksize: [int - optional] number of rows in the block rows 
    %       (necessary for the block Macaulay matrix - default = 1).
    %
    %   See also GAP, STDMONOMIALS.
    
    % MacaulayLab (2022) - Christof Vermeersch. 

    % Process the optional parameter:
    if ~exist('blocksize','var') 
        blocksize = 1;
    end
    
    % Determine gap degree and number of affine solutions:
    ma = 0;
    for k = 1:dend
        nm = blocksize*nbmonomials(k,n);
        r = sum(c<=nm);
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