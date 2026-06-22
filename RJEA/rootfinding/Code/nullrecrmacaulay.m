function [Z] = nullrecrmacaulay(Y,N,tol)
    % NULLRECRMACAULAY updates a numerical basis matrix of the null
    % space of a Macaulay matrix.
    %
    %   [Z] = NULLRECRMACAULAY(Y,N) builds a numerical basis of the 
    %   null space  of the Macaulay  matrix of degree d, using a 
    %   numerical basis of the null space Y of a lower degree Macaulay 
    %   matrix and the difference N between both Macaulay matrices.
    %
    %   Input/output variables:
    %
    %       Z: [double(q,m)] numerical basis of the updated null space.
    %       Y: [double(k,l)] numerical basis of the original null space.
    %       N: [double(p,q)] difference between both Macaulay matrices.
    %
    %   See also NULL, MACAULAYUPDATE, NULLRECRBLOCKMACAULAY.
    
    % MacaulayLab (2022) - Christof Vermeersch.
    if ~exist('tol','var') 
        tol = 10e-10;
    end

    [q,n] = size(Y); 
    X = N(:,1:q)*Y;
    [~,s,V] = svd([X N(:,q+1:end)]);
    r = nnz(s > tol);
    V = V(:,r+1:end);
    Z = [Y*V(1:n,:); V(n+1:end,:)];
end

