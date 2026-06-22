function [M] = macaulayupdate(N,problem,d,mult)
    % MACAULAYUPDATE updates the (block) Macaulay matrix.
    %
    %   [M] = MACAULAYUPDATE(N,problem,d) constructs the (block) Macaulay
    %   of a higher degree, given the (block) Macaulay matrix at the
    %   previous degree.
    %
    %   [___] = MACAULAYUPDATE(___,mult) uses a user-defined multiplication
    %   property of the polynomial basis.
    %
    %   Input argument:
    %       N [double(k,l)]: previous degree (block) Macaulay matrix.
    %       problem [struct]: system of multivariate polynomial equations
    %       or multiparameter eigenvalue problem.
    %       d [int]: degree of the (block) Macaulay matrix.
    %       mult [function - optional]: user-defined multiplication
    %       property of the polynomial basis (default = @basismon).
    %
    %   Output argument:
    %       M [double(p,q)]: sparse (block) Macaulay matrix.
    %
    %   See also MACAULAY, BASISMON, BASISCHEB.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    % Fix the important variables:
    n = problem.n; di = problem.di; s = problem.s; nnze = problem.nnze;
    p = size(problem.coef{1},2);
    q = size(problem.coef{1},3);

    % Set the default multiplication property (and compute its dilation):
    if ~exist('mult','var') 
        mult = @basismon;
    end
    [~, c] = mult(ones(1,n),ones(1,n));
    md = length(c);
    
    % Determine the number of shifts:
    shifts = zeros(s,1); % Number of shifts per polynomial.
    newshifts = zeros(s,1);
    for k = 1:s
        shifts(k) = nbmonomials(d-di{k},n);
        newshifts(k) = shifts(k) -  nbmonomials(d-di{k}-1, n);
    end
    K = monomialsmatrix(d - min(horzcat(di{:})),n);
    
    % Compute the Macaulay data:
    M = zeros(md*sum(vertcat(nnze{:}).*newshifts),3);
    pos = 0;
    row = 0;
    for k = 1:s
        eq = problem.coef{k};
        supp = problem.supp{k};
        for l = 1:size(supp,1)
            
            [re,ce,ee] = find(squeeze(eq(l,:,:)));
            [S,c] = mult(K(shifts(k)-newshifts(k)+1:shifts(k),:), supp(l,:)); 
          
            M1 = kron(ones(md,1),reshape(row + (0:newshifts(k)-1)*p + re, [length(re)*newshifts(k) 1]));
            M2 = kron(ones(md,1),reshape((position(S)-1)*q + ce, [length(ce)*newshifts(k) 1]));
            M3 = kron(c,kron(ones(newshifts(k),1),ee));
            M(pos+1:pos+newshifts(k)*md*length(ee),:) = [M1 M2 M3];
            pos = pos +newshifts(k)*md*length(ee);
        end
        row = row + p*newshifts(k);
    end
% 
% 
% 
%     row = 1;
%     pos = 1;
%     for k = 1:s
%         p = system.P{k};
%         for l = shifts(k)-newshifts(k)+1:shifts(k)
%             % Shift operation per polynomial:
%             [S,c] = mult(p(:,2:end), K(l,:)); 
%             M(pos:pos+md*nnze(k)-1,:) = [row*ones(md*nnze(k),1)  ... 
%                 position(S)' kron(c,p(:,1))];
%             row = row + 1;
%             pos = pos + md*nnze(k);
%         end
%     end
    
    % Construct the Macaulay matrix:
    orows = size(N,1);
    ocols = size(N,2);
    nrows = sum(newshifts)*p;
    ncols = nbmonomials(d,n)*q;
    M = sparse(M(:,1),M(:,2),M(:,3),nrows,ncols);
    M = [N zeros(orows, ncols-ocols); M];
end