function [Sg] = shiftmatrix(d,n,rows,shift,blocksize)
    % SHIFTMATRIX creates a shift matrix.
    %
    %   [Sg] = SHIFTMATRIX(d,n,rows,shift) creates a shift matrix in which
    %   rows are shifted by a certain shift polynomial.
    %
    %   [___] = SHIFTMATRIX(___,blocksize) considers the block version 
    %   instead.
    %
    %   Input arguments:
    %       d [int]: degree of the multivariate (block) Vandermonde matrix.
    %       n [int]: number of variables/eigenvalues.
    %       rows [int(ma)]: rows to shift.
    %       shift [double(k,n+1)]: shift polynomial.
    %       blocksize [int - optional]: size of the eigenvector.
    %
    %   Output argument:
    %       Sg [double(ma,q)]: shift matrix.
    %
    %   See also BASISMON, MONOMIALSMATRIX.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    if ~exist('blocksize', 'var')
        blocksize = 1;
    end
    k = length(rows);
    K = monomialsmatrix(d,n);
    blockpos = 1:blocksize;
    Kext = [kron(K,ones(blocksize,1)) kron(ones(size(K,1),1),blockpos')];
    l = size(Kext,1);

    Sg = zeros(k,l);
    for ki = 1:k
        pos = position(basismon(Kext(rows(ki),1:end-1),shift(:,2:end)));
        Sg(ki,(pos-1)*blocksize+Kext(rows(ki),end)) = shift(:,1);
    end
end