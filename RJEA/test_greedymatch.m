%% test_greedymatch.m
% Minimalist test to compare our greedy matching algorithm to Matlab's
% built-in matchpairs (appendix in pdf report)

addpath 'RJEA_He_Plestenjak'

% user.example = 1.5; 
if exist('user', 'var')
    S = config(user);
else
    S = config();
end

A1 = S.A1;
B1 = S.B1;
A = {A1, B1};
n = size(A{1},1); 
d = length(A); 
class_t = superiorfloat(A{:});

noise_levels = [0, 1e-12, 1e-10, 1e-8];

k = 500;
same_match = zeros(k,length(noise_levels));
for level = 1:length(noise_levels)
    noise = noise_levels(level);
    for trial = 1:k
        Xs = cell(2, 1);
        Ys = cell(2, 1);
        mus = zeros(d, 2);
        % draw multiple linear combination of A
        for iter_lc = 1:2
            % construct a normalized Gaussian vector mu
            mu = randn(d,1,class_t) + 1i*randn(d,1,class_t);
            
            mu = mu/norm(mu);
        
            % construct a linear combination A(mu)
            mus(:,iter_lc) = mu;
            M = mu(1)*A{1};
            for j = 2:d
                M = M + mu(j)*A{j};
            end
        
            % add random Gaussian noise with Frobenius norm noise to M = A(mu)
            if noise>0
                P = randn(n,class_t) + 1i*randn(n,class_t);
                P = P/norm(P,'fro');
                M = M + noise/norm(A{1})*P;
            end
        
            [X,D,Y] = eig(M);  % compute eigenvectors of M = A(mu)
        
            s = 1./diag(Y'*X); % abs(s) are condition numbers of individual eigenvalues
            S = diag(s);       % factors for 2-sided RQ, as left and right eigenvectors are normalized
            Y = Y*conj(S);      % ensures Y'*X = Id
            Xs{iter_lc} = X;
            Ys{iter_lc} = Y;
        
        end
        
        % combine eigenvectors found with matching function
        p_greedy = greedyMatching(Xs);
        S = abs(Xs{1}'*Xs{2});
        M = matchpairs(-S, 1e10);   
        p_hungarian = zeros(n, 1);
        p_hungarian(M(:,1)) = M(:,2);
        same_match(trial, level) = all(p_greedy{1} == p_hungarian);
    end
    fprintf('For noise level %5.2e, greedy agrees with Hungarian: %.2f%%\n', noise_levels(level), 100*mean(same_match(:, level)));
end

