function [D, commuting_family] = shiftnullspace(Z,G,dgap,selectedrows,blocksize,mult)
    % SHIFTNULLSPACE solves a system of multivariate polynomial equations
    % or multiparameter eigenvalue problem by shifting the null space of a
    % corresponding (block) Macaulay matrix. The shift problem is solved 
    % via multiple Schur decompositions.
    %
    %   [D] = SHIFTNULLSPACE(Z,G,dgap,selectedrows) solves a system of
    %   multivariate polynomial equations via basis matrix of the null 
    %   space of a corresponding Macaulay matrix.
    %
    %   [D] = SHIFTNULLSPACE(___,blocksize) solves a multiparameter
    %   eigenvalue problem via a basis matrix of the null space of a 
    %   corresponding block Macaulay matrix.
    %
    %   [D] = SHIFTNULLSPACE(___,blocksize,mult) uses a specific polynomial     
    %   basis multiplication property to perform the shifts.
    %
    %   Input/output variables:
    %
    %       D: [cell(n+1)] evaluation of the shift polynomial in the affine
    %       solutions + the affine solutions (last n elements of D).
    %       Z: [double(q,ma)] basis matrix of the corresponding (block) 
    %       Macaulay matrix.
    %       G: [int(l,n+1)] shift polynomial of the first shift problem.
    %       dgap: [int] degree of first gap block.
    %       selectedrows: [int(ma,1)] indices of affine standard monomials 
    %       in the basis matrix - if set to NaN, then a blocked approach is 
    %       used.
    %       blocksize: [int - optional] the size of a block row in the 
    %       basis matrix (1 corresponds to a system of multivariate 
    %       polynomial equations, while > 1 corresponds to a multiparameter
    %       eigenvalue problem - default = 1).
    %       mult: [function] polynomial basis multiplication property
    %       (default = @basisvander).
    %
    %   See also MONOMIALSMATRIX, SCHUR.
    
    % MacaulayLab (2022) - Christof Vermeersch.
    
    % Set the default block size (SMPE):
    if ~exist('blocksize','var') 
        blocksize = 1;
    end

    % Set the default multiplication property (and compute its dilation):
    if ~exist('mult','var') 
        mult = @basismon;
    end
    n = size(G,2) - 1;
    [~, c] = mult(ones(1,n),ones(1,n));
    md = length(c);    
    
    % Create large-enough matrix with monomials:
    K = monomialsmatrix(dgap - 1,n);

    if isnan(selectedrows) % Block-wise approach.
        % Determine the selected rows to shift:
        selectedrows = blocksize*nbmonomials(dgap - 1,n);
        SZ = cell(n+2,1);
        
        % Perform the row selections:
        SZ{1} = Z(1:selectedrows,:);

        % Perform the row combinations:
        for s = 1:n+1
            % Construct the shift polynomial:
            if s == 1
                H = G;
            else
                H = [1 zeros(1,s-2) 1 zeros(1,n-s+1)];
            end

            % Construct corresponding shift matrix:
            S = zeros(selectedrows,size(Z,2));
            for k = 1:size(H,1)
                [D,c] = mult(K,H(k,2:end));
                d = position(D);
                blockvector = 1:blocksize;
                for l = 1:md
                    shiftrows = (l-1)*nbmonomials(dgap - 1,n)+1:l*nbmonomials(dgap - 1,n);
                    hitrows = (d(shiftrows)-1)*blocksize + blockvector';
                    S = S + c(l)*H(k,1)*Z(hitrows(:),:);
                end
            end
            SZ{s+1} = S;
        end
    else % Row-wise approach.
        SZ = cell(n+2,1);
        % Perform the row selections:
        SZ{1} = Z(selectedrows,:);
        
        % Perform the row combinations:
        for s = 1:n+1
            % Construct the shift polynomial:
            if s == 1
                H = G;
            else
                H = [1 zeros(1,s-2) 1 zeros(1,n-s+1)];
            end
            
            % Construct corresponding shift matrix:
            S = zeros(length(selectedrows),size(Z,2));
            
            % Construct corresponding shift matrix:
            blockvector = 1:blocksize;
            Kextended = kron(K,ones(blocksize,1));
            Lextended = kron(ones(size(K,1),1),blockvector');
            for k = 1:size(H,1)
                [D,c] = mult(Kextended(selectedrows,:),H(k,2:end));
                d = position(D);
                for l = 1:md
                    shiftblocks = (l-1)*length(selectedrows)+1:l*length(selectedrows);
                    hitblocks = d(shiftblocks)';
                    hitrows = (hitblocks - 1)*blocksize + Lextended(selectedrows);
                    S = S + c(l)*H(k,1)*Z(hitrows(:),:);
                end
            end
            SZ{s+1} = S;

        end
    end
    
    % Solve the first eigenvalue problem:
    A = SZ{1};
    B = SZ{2};
    [Q,S] = schur(pinv(A)*B,'complex');
    D = cell(n+1,1);
    D{1} = diag(S);
    commuting_family=cell(n+1,1);
    commuting_family{1} = pinv(A)*B;
    % Solve the other eigenvalue problems:
    for s = 2:n+1
        B = SZ{s+1};
        commuting_family{s} = pinv(A)*B;
        S = Q'*pinv(A)*B*Q;
        D{s} = diag(S);
    end
end