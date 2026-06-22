function [Z] = columncompression(Z,dcompr,n,ma,blocksize)
    % COLUMNCOMPRESSION compresses a basis matrix of the (right) null space
    % of a (block) Macaulay matrix.
    %
    %   [Z] = COLUMNCOMPRESSION(Z,dcompr,n,ma) compresses a basis matrix of
    %   the (right) null space of a (block) Macaulay matrix with
    %   block size equal to 1.
    %
    %   [___] = COLUMNCOMPRESSION(___,blocksize) considers a block size 
    %   equal to blocksize.
    %
    %   Input arguments:
    %       Z [double(k,mb)]: basis matrix of the (right) null space of a
    %       (block) Macaulay matrix.
    %       dcompr [int]: degree of the compression.
    %       n [int]: number of variables.
    %       ma [int]: number of columns after performing the compression.
    %       blocksize [int - optional]: size of a block (default = 1).
    %
    %   Output argument:
    %	    Z [double(l,ma)]: basis matrix of the compressed (right) null 
    %       space of a (block) Macaulay matrix.

    % MacaulayLab (2023) - Christof Vermeersch.

    % TODO: use more simple column compression.
    
    % Set the default block size:
    if ~exist('blocksize','var') 
        blocksize = 1;
    end

    nbrows = blocksize*nbmonomials(dcompr,n);
    [~,~,Q] = svd(Z(1:nbrows,:));
    V = Z*Q;
    Z = V(1:nbrows,1:ma);
end