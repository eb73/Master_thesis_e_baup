function [S] = config(user)
% CONFIG(user) initializes all variables needed to run RJEA (matrices,
% number of trials, method, auxiliary parameters)
%
% Input:
%   user: struct containing user specified settings (to replace default
% values)
%
% Output:
%   S: struct containing all parameters needed for the main (to call RJEA & variants,
% as well as some parameters for plot settings)

%% general setting
if nargin == 0 || isempty(user)
    user = struct();
end

% default values
S.example = 1.1;
S.method = 'best_res'; % other options: 'ref', 'mean', 'smaller_subspace'
S.Nsamples = 1000;
S.download_output = false;
S.nb_trials = 2;
S.rq_type = 'one-sided'; % other option: 'two-sided'


% override with user values
fn = fieldnames(user);
for k = 1:numel(fn)
    S.(fn{k}) = user.(fn{k});
end

% in case use different name for rq1/rq2
if strcmp(S.rq_type, 'rq1')
    S.rq_type = 'one-sided';
elseif strcmp(S.rq_type, 'rq2')
    S.rq_type = 'two-sided';
end

if ismember(S.method, {'best residual', 'best_residual'})
    S.method = 'best_res';
elseif ismember(S.method, {'smaller subspace', 'subspace'})
    S.method = 'smaller_subspace';
end

%% example computation

uniform_conditioned = 1;
class_t = 'double';

example = S.example;
if (example > 0) && (example < 3)
    if (example > 0) && (example < 1)
        % example 5.1 from RJEA, Plestenjak, He, Kressner
        da = [1 1 1 2 2 2 3]';  
        db = [1 2 3 1 2 3 3]';
        DA = diag(da);
        DB = diag(db);
        noise_level = [-inf -14 -12 -10]; % [-14, -10];%  
        kappa = 1e2;
        index = [1];
        string_eig = '(1,1)';
        if example == 0.1
            plot_figures = [1, 2];
        elseif example == 0.2
            noise_level = [-12 -11 -10 -9];
            plot_figures = [10];    
        end
        ylim_err = [1e-15, 1e-8];
        ylim_res = [-18, -10];
    elseif (example > 2) && (example < 3)
        % example 5.3 from RJEA, Plestenjak, He, Kressner
        da = [1 1 1 2 2 2 3]';  
        db = [1 2 3 1 2 3 3]';
        DA = diag(da);
        DB = diag(db);
        noise_level = [-inf -12 -10 -8];
        kappa = 1e4;  % condition number of eigenvalue matrix X (try 1e2, 1e4, 1e6, 1e8)
        uniform_conditioned = 0;
       plot_figures = [1, 2]; 
        if example == 2.1
            index = 1; string_eig = '(1,1)';
        elseif example == 2.2
            index = 4; string_eig = '(2,1)';
        end
        ylim_err = [1e-14, 1e-4];
        ylim_res = [-18, -10];
    elseif (example > 1) && (example < 2)
        n = 10; % change if needed, requires multiple of 5
        if example == 1.2
            n = 50;
        elseif example == 1.3
            n = 100;
        end
        
        da = repelem((1:n/5), 5)';
        db = repmat((1:5),1,n/5)'; %(1:n)'; 
        DA = diag(da);
        DB = diag(db);
        ylim_err = [1e-16, 1e-6];
        ylim_res = [-18, -9];
        noise_level = [-inf -12 -10 -8];
        if example == 1.5
            eps_ex = 1e-10;
            db(2) = 1-eps_ex;
            db(3) = 1+eps_ex;
            %db = repmat([1, 1+1e-14, 1+1e-12, 1+1e-10, 1+1e-8], 1, 2)'; %(1:n)'; 
            DB = diag(db);
            ylim_err = [1e-15, 1e-4];
            ylim_res = [-18, -8];
            noise_level = [-Inf, -14];%[-14 -12 -10 -8 -6];
        end
        kappa = 1e2;
        index = [1];
        string_eig = '(1,1)';
        plot_figures = [1, 2];
    end
    n = length(da);
    rng(1)
    if kappa == 1
        X1 = eye(n,class_t);
    else
        if uniform_conditioned == 1
            X1 = condmatX(n,kappa,class_t); 
        else
            % first part with two vectors is well conditioned
            X2 = condmatX(n-2,kappa,class_t); 
            X1 = zeros(n,class_t); X1(1:2,1:2) = eye(2,class_t); X1(3:n,3:n) = X2; % ill conditioned X but small kappa1
            X1 = randn(n,class_t)*X1;
        end
    end
    Y1 = inv(X1);
    A1 = X1*DA*Y1;
    B1 = X1*DB*Y1;
    
elseif example == 4.1
    plot_figures = [1, 2];
    index = 1;
    string_eig = '(3+i,12)';
    A1 = ones(3) + diag(1i*ones(3,1));
    B1 = [7 0 5; 2 4 6; 3 8 1];
    da = [3+1i; 1i; 1i]; % beware do not use ' here because will conjugate i
    db = [12 -2*sqrt(6) 2*sqrt(6)]';
    n = length(da);
    noise_level = [-inf -12];
    %true_je = [1 5 5; 1 2*(1+sqrt(6)) 2*(1-sqrt(6)); 1 -(7+2*sqrt(6)) -(7-2*sqrt(6))];
    ylim_err = [-Inf, Inf];
    ylim_res = [-Inf, Inf];
elseif example == 4.2
    plot_figures = [1, 2];
    index = 1;
    string_eig = '(2, -1)';
    A1 = [0 4 0 0; 
          1 0 0 0;
          0 0 0 4;
          0 0 1 0];
    B1 = [0 0 1 0;
          0 0 0 1;
          1 0 0 0;
          0 1 0 0];
    da = [2 2 -2 -2]';
    db = [-1 1 -1 1]';
    n = length(da);
    noise_level = [-inf -12 -10 -8];
    %true_je = [2 2 2 2; 1 1 -1 -1; -2 2 -2 2; -1 1 1 -1];
    ylim_err = [1e-18, 1e-7];
    ylim_res = [-18, -8];

end

Aset = {A1,B1};
exact_eig = [da db];
noise = 10.^noise_level;
noise_for_Amu = 0;

center_eig = sum(exact_eig(index,:),1)/length(index);
if length(index)==1
    delta = 0;
else
    delta = norm(exact_eig(index(1),:)-center_eig);
    for j = 2:length(index)
        delta = max(delta, norm(exact_eig(index(j),:)-center_eig));
    end
end
delta = delta*sqrt(length(index));

d = length(Aset);
realfamily = 1;
for j = 1:d
    realfamily = and(realfamily,isreal(Aset{j}));
end

gapA = Inf;
for j = 1:n
    if j~=index
        gapA = min(norm(exact_eig(j,:)-exact_eig(index,:)),gapA);
    end
end

min_noise = norm([norm(A1) norm(B1)])*eps(class_t)/2;

% random noise matrices we add to A to simmulate near commutativity
PA = cell(d,1);
for j = 1:d 
    if realfamily
        PA{j} = randn(n,class_t);
    else
        PA{j} = randn(n,class_t) + 1i*randn(n,class_t);
    end
    PA{j} = PA{j}/norm(PA{j},'fro');
end

has_subspace_rq2 = strcmp(S.method, 'smaller_subspace') && ismember(S.rq_type, {'two-sided', 'both'});
%has_mean_ex1_5 = strcmp(S.method, 'mean') && example == 1.5;

if has_subspace_rq2 && S.nb_trials == 2
    ylim_err = [ylim_err(1), Inf];
    ylim_res = [ylim_res(1), Inf];
end

%% assign all computed values to structure

S.A1 = A1;
S.B1 = B1;
S.Aset = Aset;
S.PA = PA;
S.noise = noise;
S.noise_level = noise_level;
S.noise_for_Amu = noise_for_Amu;
S.min_noise = min_noise;
S.delta = delta;
S.gapA = gapA;
S.index = index;
S.d = d;
S.exact_eig = exact_eig;

S.plot_figures = plot_figures;
S.string_eig = string_eig;
S.ylim_err = ylim_err;
S.ylim_res = ylim_res;

end