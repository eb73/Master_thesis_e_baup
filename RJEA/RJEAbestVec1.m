function [best_z, best_res] = RJEAbestVec1(B, complex, n_trials)
% Performs a modified version of the RJEA algorithm, 
% and return the best eigenvector in terms of 1S residual 
%
% Input
%  - B: set of (nearly) commuting 2x2 matrices
%  - complex: boolean, true if use complex values, false o/w
%  - n_trials: number of runs
% Output
%  -  best_z: best eigenvector
%  -  best_res: best associated residual

class_t = superiorfloat(B{:});
d = length(B);   
best_res = Inf;
best_z = eye(size(B{1},1), 1);   % first canonical basis vector

for j = 1:d
    norms(j) = norm(B{j});
end

realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(B{j}));
end

for it = 1:n_trials
% construct a normalized Gaussian vector mu
    if complex
        mu = randn(d,1,class_t) + 1i*randn(d,1,class_t);
    else
        mu = randn(d,1,class_t);
    end
    mu = mu/norm(mu);
    
    % construct a linear combination A(mu)
    M = mu(1)*B{1};
    for j = 2:d
       M = M + mu(j)*B{j};
    end
 
    [X,~,~] = eig(M);  % compute eigenvectors of M = B(mu)
    for idx = 1:size(X,2)
        z = X(:,idx);
        z = z/norm(z);
        % compute residual 
        comp_eig = zeros(1,d,class_t);
        for j = 1:d
            comp_eig(:,j) = z'*B{j}*z / norm(z);
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
        end
    end

end
