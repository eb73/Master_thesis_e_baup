function [M] = macaulay(problem,d,mult)
    % MACAULAY constructs the (block) Macaulay matrix.
    %
    %   [M] = MACAULAY(problem,d) constructs the (block) Macaulay matrix of
    %   a given degree in the standard monomial basis.
    %
    %   [___] = MACAULAY(___,mult) uses a user-defined multiplication
    %   property of the polynomial basis.
    %
    %   Input argument:
    %       problem [struct]: system of multivariate polynomial equations
    %       or multiparameter eigenvalue problem.
    %       d [int]: degree of the (block) Macaulay matrix.
    %       mult [function - optional]: user-defined multiplication
    %       property of the polynomial basis (default = @basismon).
    %
    %   Output argument:
    %       M [double(p,q)]: sparse (block) Macaulay matrix.
    %
    %   See also BASISMON, BASISCHEB.
    
    % MacaulayLab (2023) - Christof Vermeersch.
    
    % Fix the important variables:
    [s,n,di,~,nnze,p,q] = probleminfo(problem);

    % Set the default multiplication property (and compute its dilation):
    if ~exist('mult','var') 
        mult = @basismon;
    end
    [~, c] = mult(ones(1,n),ones(1,n));
    md = length(c);

    % Determine the number of shifts (vectorized not possible):
    shifts = zeros(s,1); % Number of shifts per equation.
    for k = 1:s
        shifts(k) = nbmonomials(d - di{k},n);
    end
    
    % Create the monomial vectors:
    K = monomialsmatrix(d - min(horzcat(di{:})),n);
    
    % Compute the Macaulay data:
    nnzm = vertcat(nnze{:}).*shifts;
    M = zeros(md*sum(nnzm),3);
    pos = 0;
    row = 0;
    for k = 1:s
        eq = problem.coef{k};
        supp = problem.supp{k};
        for l = 1:size(supp,1)
            
            [re,ce,ee] = find(squeeze(eq(l,:,:)));
            [S,c] = mult(K(1:shifts(k),:), supp(l,:)); 

            M1 = kron(ones(md,1),reshape(row + (0:shifts(k)-1)*p + re, [length(re)*shifts(k) 1]));
            M2 = reshape((position(S)-1)*q + ce, [md*length(ce)*shifts(k) 1]);
            M3 = kron(c,kron(ones(shifts(k),1),ee));
            M(pos+1:pos+shifts(k)*md*length(ee),:) = [M1 M2 M3];
            pos = pos +shifts(k)*md*length(ee);
        end
        row = row + p*shifts(k);
    end
    
    % Construct the Macaulay matrix:
    nrows = sum(shifts)*p;
    ncols = nbmonomials(d,n)*q;
    M = sparse(M(:,1),M(:,2),M(:,3),nrows,ncols);
end