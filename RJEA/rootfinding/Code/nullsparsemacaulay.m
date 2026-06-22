function [Z] = nullsparsemacaulay(Z,problem,d,tol,mult)
    % NULLSPARSEMACAULAY computes a numerical basis matrix of the null
    % space without building a Macaulay matrix.
    %
    %   [Z] = NULLSPARSEMACAULAY(Z,system,d) computes a numerical basis of 
    %   the null space of the Macaulay matrix of degree d that comprises 
    %   system, without building the corresponding Macaulay matrix.
    %
    %   [Z] = NULLSPARSEMACAULAY(___,mult) uses a specific polynomial basis
    %   multiplication property to perform the shifts.
    %
    %   Input/output variables:
    %
    %       Z: [double(q,mb)] numerical basis of the null space.
    %       system: [struct] system of multivariate polynomial equations.
    %       d: [int] degree of the Macaulay matrix.
    %       mult: [function] polynomial basis multiplication property
    %       (default = @basisvander).
    %
    %   See also NULL, NULLRECRROW, NULLRECRBLOCKMACAULAY, 
    %   NULLSPARSEBLOCKMACAULAY, BLOCKMACAULAYUPDATE.
    %
    %   Note: this implementation is quite naive and receives an update in 
    %   one of the next releases.
    
    % MacaulayLab (2022) - Christof Vermeersch. 
    
    % Retrieve the information of the SMPE:
    [s,n,di,~,~,p,q] = probleminfo(problem);

    if ~exist('tol','var') 
        tol = 10e-10;
    end

    % Set the default multiplication property (and compute its dilation):
    if ~exist('mult','var') 
        mult = @basismon;
    end
    [~, c] = mult(ones(1,n),ones(1,n));

    % Determine the number of shifts:
    newshifts = zeros(s,1); % Number of new shifts per polynomial..
    for k = 1:s
        newshifts(k) = nbmonomials(d - di{k},n) ...
            - nbmonomials(d - di{k} - 1,n);
    end
    lastcol = nbmonomials(d-1,n)*q;
    newcolumns = (nbmonomials(d,n) - nbmonomials(d-1,n))*q;
    newrows = cumsum([0; newshifts]);
    nul = size(Z,2);
    
    % Create the monomial vectors:
    K = monomialsmatrix(d,n);

    % Increase the basis matrix:
    
    X = zeros(sum(newshifts),nul);
    Y = zeros(sum(newshifts),newcolumns);
    E = eye(newcolumns);
    for k = 1:s
        eq = problem.coef{k};    % Concatenate of all the polynomials? (removes one loop) in the MEP:
        supp = problem.supp{k};
        T = kron(c',reshape(shiftdim(eq,1),[p,q*size(eq,1)])); % Elements per polynomial.
        for l = 1:newshifts(k)
            shift = nbmonomials(d - di{k} - 1,n) + l;
            S = mult(supp,K(shift,:)); 
            S = (position(S)-1)*q + (1:q)';
            idx = (S(S(:) < lastcol + 1)); % Determine existing indices.
            idxval = S(:) < lastcol + 1;
            idy = (S(S(:) > lastcol)); 
            idyval = S(:) > lastcol;
            X((l-1)*p+1+newrows(k):l*p+newrows(k),:) = T(:,idxval)*Z(idx(:),:);
            Y((l-1)*p+1+newrows(k):l*p+newrows(k),:) = T(:,idyval)*E(idy(:)-lastcol,:); 
        end
    end

    % Update the numerical basis matrix:
    [~,s,V] = svd([X Y]);
    r = nnz(s > tol);
    V = V(:,r+1:end);
    Z = [Z*V(1:nul,:); V(nul+1:end,:)];
end