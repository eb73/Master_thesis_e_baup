function [D] = shiftcolumnspace(M,G,dgap,caffine,cinfinity,blocksize,mult)
    % SHIFTCOLUMNSPACE solves a system of multivariate polynomial equations
    % by shifting the column space of a corresponding Macaulay matrix. The 
    % shift problem is solved via multiple Schur decompositions.
    %
    %   [D] = SHIFTCOLUMNSPACE(M,G,dgap,caffine,cinfinity) solves a system 
    %   of multivariate polynomial equations via the Macaulay matrix M.
    %
    %   [D] = SHIFTCOLUMNSPACE(___,mult) solves a system 
    %   of multivariate polynomial equations via the Macaulay matrix M.
    %
    %   Input/output variables:
    %
    %       D: [cell(n+1)] evaluation of the shift polynomial in the affine
    %       solutions + the affine solutions (last n elements of D).
    %       M: [double(p,q)] Macaulay matrix
    %       G: [int(l,n+1)] shift polynomial of the first shift problem.
    %       dgap: [int] degree of first gap block.
    %       caffine: [int(ma,1)] indices of affine standard monomials 
    %       cinfinity: [int(minfinity,1) indices of the standard monomials
    %       related to solutions at infinity.
    %       mult: [function] polynomial basis multiplication property
    %       (default = @basisvander).
    %
    %   See also MONOMIALSMATRIX, SCHUR.
    
    % MacaulayLab (2022) - Christof Vermeersch.

%     warning('The column space based shift is not optimized for efficiency')
    if ~exist('blocksize','var') 
        blocksize = 1;
    end

    % Set the default multiplication property (and compute its dilation):
    if ~exist('mult','var') 
        mult = @basismon;
    end
    [~, c] = mult(1,1);
    md = length(c);

    % Determine the number of variables and affine solutions:
    n = size(G,2) - 1;
    ma = length(caffine);
    
    % Create large-enough matrix with monomials:
    K = monomialsmatrix(dgap - 1,n);
    blockvector = 1:blocksize;
    Kextended = kron(K,ones(blocksize,1));
    Lextended = kron(ones(size(K,1),1),blockvector');
    SM = cell(n+1,1);
    
    % Perform the shifts:
    for s = 1:n+1
        % Construct the shift polynomial:
        if s == 1
            H = G;
        else
            H = [1 zeros(1,s-2) 1 zeros(1,n-s+1)];
        end
        
        % Construct corresponding shift matrix:
        cshift = zeros(length(caffine),md*size(H,1));
        for k = 1:ma
            [D,~] = mult(Kextended(caffine(k),:),H(:,2:end));
            blockposition = (position(D)-1)*blocksize;
            cshift(k,:) = blockposition+Lextended(caffine(k));
        end
        chit = setdiff(cshift,caffine);
        cremainder = setdiff(setdiff(setdiff(1:size(M,2),caffine), ...
            cshift),cinfinity);

        shiftvalue = kron(c,H(:,1));

        % Construct corresponding selection matrix:
        Sg = zeros(ma,ma+length(chit));
        for k = 1:ma
            [~,pos,posinshift] = intersect(caffine,cshift(k,:));
            if size(pos,1) ~= 0
                Sg(k,pos) = shiftvalue(posinshift);
            end
            [~,pos,posinshift] = intersect(chit,cshift(k,:));
            pos = pos + ma;
            if size(pos,1) ~= 0
                Sg(k,pos) = shiftvalue(posinshift);
            end
        end
        
        % Construct R33 and R34:
        M1 = M(:,caffine);
        M2 = M(:,chit);
        M3 = M(:,cremainder);
        
        [~,R] = qr(flip([M1 M2 M3],2),0);
        R = flip(R,2);
        R33 = R(size(M3,2)+1:size(M3,2)+size(M2,2),size(M1,2)+ ...
            1:size(M1,2)+size(M2,2));
        R34 = R(size(M3,2)+1:size(M3,2)+size(M2,2),1:size(M1,2));
        
        % Build second coefficient matrix of the eigenvalue problems:
        SM{s} = Sg*[eye(ma);-inv(R33)*R34];
    end

    % Solve the first eigenvalue problem:
    A = SM{1};
    
    [Q,S] = schur(A,'complex'); 
    D = cell(n+1,1);
    D{1} = diag(S);
      
    % Solve the other eigenvalue problems:
    for s = 2:n+1
        A = SM{s};
        
        S = Q'*A*Q;
        D{s} = diag(S);
    end
end
