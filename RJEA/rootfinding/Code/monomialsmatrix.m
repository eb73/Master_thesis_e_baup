function [K] = monomialsmatrix(d,n,order)
    % MONOMIALSMATRIX creates a matrix with the monomial vectors.
    %
    %   [K] = MONOMIALSMATRIX(d,n) creates a matrix with the monomial 
    %   vectors of the monomial basis using the graded inverse 
    %   lexicographic monomial ordering.
    %
    %   [___] = MONOMIALSMATRIX(___,order) uses a user-defined monomial
    %   ordering.
    %
    %   Input arguments:
    %       d [int]: maximum total degree of the monomials.
    %       n [int]: number of variables of the monomials.
    %       order [function - optional]: user-defined monomial ordering
    %       (default = @posgrinvlex).
    %
    %   Output argument: 
    %       K: [double(d,n)] monomial basis in its vector representation.
    %
    %   See also NBMONOMIALS, POSITION.
        
    % MacaulayLab (2023) - Christof Vermeersch.
    
    if d == 0
        K = zeros(1,n);
    else
        I = eye(n,n);
        K = zeros(nbmonomials(d,n),n);
        K(2:n+1,:) = I;
        blocklength = ones(n,1);
        endpoint = 1;
        for di = 2:d
            blocklength = cumsum(blocklength, 'reverse');
            blockpos = [0; cumsum(blocklength)];
            endpoint = endpoint + blocklength(1);
            for ni = 1:n
                K(endpoint+blockpos(ni)+1:endpoint+blockpos(ni+1),:) = ...
                    K(endpoint-blocklength(ni)+1:endpoint,:) + I(ni,:);             
            end
        end
    end

    % Order the monomials according to a specific monomial ordering.
    if exist('order','var') 
        K = K(order(K),:);
    end
end