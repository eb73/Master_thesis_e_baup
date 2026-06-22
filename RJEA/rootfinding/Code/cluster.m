function [D,centers,elements] = cluster(D,tol)
    % CLUSTER determines the different clusters of solutions candidates in 
    % order to impove their accuracy.
    %
    %   [D] = CLUSTER(D) clusters the different solution candidates in D
    %   based on the first component in D.
    %
    %   [___,centers,elements] = CLUSTER(D) gives the different centers
    %   and elements of every identified cluster.
    %
    %   [___] = CLUSTER(___,tol) uses a user-specified clustering
    %   tolerance.
    %
    %   Input arguments:
    %       D [cell(n+1,1)]: different components of the affine solution
    %	    candidates before clustering - the clustering algorithm
    %	    uses the first element of D to determine the clusters.
    %       tol [double - optional]: the tolerance to determine when two 
    %	    elements belong to the same cluster (default = 1e-10).
    %
    %   Output arguments:
    %	    D [cell(n+1,1)]: different components of the affine solution
    %	    candidates after clustering.
    %	    centers [double(nc,1)]: cluster centers. 
    % 	    elements [cell(nc,1)]: the indices of the solution candidates
    % 	    that belong to an identified cluster.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    % TODO: allow other data structures.
    
    % Process the optional argument:
    if ~exist('tol','var') 
        tol = 1e-10;
    end
    
    % Cluster the solutions based on the values of the first element of D:
    S = D{1};
    n = length(D) - 1;
    N = length(S);
    centers = zeros(N,1);
    elements = cell(N,1);
    centers(1) = S(1);
    elements{1} = 1;
    nc = 1;
    for k = 2:N
        dist = zeros(nc,1);
        for l = 1:nc
            centers(l) = mean(S(elements{l}));
            dist(l) = abs(S(k)-centers(l));
        end
        [m,ind] = min(dist);
        if m < tol
            elements{ind} = [elements{ind}; k];
        else
            nc = nc + 1;
            centers(nc) = S(k);
            elements{nc} = k;
        end
    end
    centers = centers(1:nc);
    elements = elements(1:nc);
    
    % Re-order the components of the solutions based on the clustering 
    % results:
    for k = 0:n
        S = D{k+1};
        Sclustered = zeros(length(S),1);
        ca = 0;
        for ci = 1:length(centers)
            E = S(elements{ci});
            ca = ca + length(E);
            Sclustered(ca-length(E)+1:ca) = mean(E)*ones(length(E),1);
        end
        D{k+1} = Sclustered;
    end
end     