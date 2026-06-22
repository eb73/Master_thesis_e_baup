function [Z] = nullrecrrow(Y,N,tol)
    % NULLRECROW computes a numerical basis matrix of the null space of a
    % block row matrix.
    %
    %   [Z] = NULLRECRROW(Y,N) builds a numerical basis Z of the null space 
    %   of a (block) row  matrix of degree d, using a numerical basis Y of 
    %   the null space of a lower degree (block) row matrix and the 
    %   difference N between both (block) row matrices.
    %
    %   [Z] = NULLRECRROW(___,tol) uses a user-specified tolerance to catch
    %   near-zero matrices.
    %
    %   Input/output variables:
    %
    %       Z: [double(q,m)] numerical basis of the updated null space.
    %       Y: [double(k,l)] numerical basis of the original null space.
    %       N: [double(p,q)] difference between both (block) row matrices.
    %
    %   See also NULL.
    
    % MacaulayLab (2022) - Christof Vermeersch. 
    
    % Process the optional parameter:
    if ~exist('tol','var') 
        tol = 1e-10;
    end

    X = N*Y;
    [~,S,V] = svd(X);
    r = nnz(diag(S) > tol);
    if r == 0
        V = eye(size(V,1));
    end
    Z = Y*V(:,r+1:end);
end