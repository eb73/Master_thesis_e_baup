function [best_z, best_w, best_res] = RJEAbestVec2(B, complex, n_trials)
% Performs a modified version of the RJEA algorithm, 
% and returns the best eigenvectors in terms of 2S residual 
%
% Input
%  - B: set of (nearly) commuting 2x2 matrices
%  - complex: boolean, true if use complex values, false o/w
%  - n_trials: number of runs
% Output
%  -  best_z: best right eigenvector
%  -  best_w: best left eigenvector
%  -  best_res: best associated residual

class_t = superiorfloat(B{:});
d = length(B);    % number of commuting matrices in set A 
% initialize best residual, right/left eigenvec
best_res = Inf;
best_z = eye(size(B{1},1), 1);   % first canonical basis vector
best_w = eye(size(B{1},1), 1);
for j = 1:d
    norms(j) = norm(B{j});
end
realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(B{j}));
end

for it = 1:n_trials
    if complex
        mu = randn(d,1,class_t) + 1i*randn(d,1,class_t);
    else
        mu = randn(d,1,class_t);
    end
    mu = mu/norm(mu);
    
    % construct a linear combination B(mu)
    M = mu(1)*B{1};
    for j = 2:d
       M = M + mu(j)*B{j};
    end
 
    [X,~,Y] = eig(M);  % compute eigenvectors of M = B(mu)
    
    %lamu = diag(D); % eigenvalue of B(mu)
    s = 1./diag(Y'*X); % abs(s) are condition numbers of individual eigenvalues
    S = diag(s);       % factors for 2-sided RQ, as left and right eigenvectors are normalized   
    Y1 = Y*conj(S);    % left eigenvectors scaled so that Y1'*X = Id
    
    for idx = 1:size(X,2)
        z = X(:,idx);
        z = z/norm(z);% shouldn't be necessary
        w = Y1(:, idx); % take associated left eigenvector
        if abs(w'*z - 1) > 1e-10
            warning('problem')
        end
        % compute residual 
        comp_eig = zeros(1,d,class_t);
        for j = 1:d
            comp_eig(:,j) = (w' * B{j} * z) / (w' * z);
        end
        tmpres = compute_residual(B,comp_eig,norms);
        r2 = tmpres(:,1).^2;
        for j = 2:d
            r2 = r2 + tmpres(:,j).^2;
        end
        res = norm(sqrt(r2));   

        if res < best_res
            best_res = res;
            best_z = z;
            best_w = w;
        end
    end

end