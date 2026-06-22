function [V] = vandermonde(d,n,x,z)
    % VANDERMONDE constructs the multivariate (block) Vandermonde basis
    % matrix.
    %
    %   [V] = VANDERMONDE(d,n,x) constructs the multivariate Vandermonde
    %   basis matrix in the standard monomial basis.
    %
    %   [___] = VANDERMONDE(___,z) considers the block version instead.
    %
    %   Input arguments:
    %       d [int]: degree of the multivariate (block) Vandermonde matrix.
    %       n [int]: number of variables/eigenvalues.
    %       x [double(ma,n)]: solutions of the problem.
    %       z [double(ma,blocksize) - optional]: eigenvectors of the
    %       multiparameter eigenvalue problem.
    %
    %   Output argument:
    %       V [double(q,ma)]: multivariate (block) Vandermonde basis
    %       matrix.
    %
    %   See also BASISMON, MONOMIALSMATRIX.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    if ~exist('z','var')
        z = ones(size(x,1),1);
    end
    
    K = monomialsmatrix(d,n);
    V = zeros(size(K,1)*size(z,2),size(x,1));
    for k = 1:size(x,1)
        V(:,k) = kron(prod(x(k,:).^K,2),z(k,:).');
    end
end