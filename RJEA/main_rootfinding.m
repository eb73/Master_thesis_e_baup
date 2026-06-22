%% main_rootfinding.m
% Main file used for the rootfinding example (RJE application)

addpath(genpath('rootfinding'))

%% Load problem

if ~exist('example', 'var')
    example = 'redeco8';
end


if strcmp(example, 'redeco8')
    problem = redeco8();
    load("redeco8.mat",'A')
elseif strcmp(example, 'rose')
    load('rose.mat', 'A')
    problem = rose();
elseif strcmp(example, 'noon5')
    load("noon5.mat")
    problem = noon5();
elseif strcmp(example, 'random')
    warning('This example may take a long time to run')
    problem = randomsystem(3,10,3);
    options = struct();
    options.recursive = 'recursive';
    [commonroots,options,A] = macaulaylab(problem, 30, options);
else
     error('Example not supported')
end

repeats = 1;
trial_values = 2:6;
nTrials = length(trial_values);

fprintf('Running rootfinding with example %s, %d repeats and max k=%d \n', example, repeats, nTrials+1)

%% Run algorithms

% Reference (rjd)
err_rjd = zeros(repeats,1);
time_rjd = 0;

for r = 1:repeats

    fprintf("RJD repeat #%d\n", r);

    start = tic;
    rjd_sol = rjd(A);
    time_rjd = time_rjd + toc(start);

    [err_rjd(r),~] = residuals(problem,rjd_sol,@evalmon);

end

time_rjd = time_rjd/repeats;

% My methods (run nTrials time, each time with different k)
err_mean1   = zeros(nTrials,repeats);
err_mean2   = zeros(nTrials,repeats);
err_bestres = zeros(nTrials,repeats);
err_sub1    = zeros(nTrials,repeats);
err_sub2    = zeros(nTrials,repeats);

time_mean    = zeros(nTrials,1);
time_bestres = zeros(nTrials,1);
time_sub1    = zeros(nTrials,1);
time_sub2    = zeros(nTrials,1);

d = length(A);
normsA = zeros(d);
for j = 1:d
    normsA(j) = norm(A{j});
end

for k = 1:nTrials

    nb_trials = trial_values(k);
    fprintf("Nb trial k=%d \n", nb_trials)

    tmean = 0;
    tbest = 0;
    tsub1 = 0;
    tsub2 = 0;

    for r = 1:repeats
        fprintf("Repeat #%d\n", r);

        % Mean
        start = tic;
        [rjdmean_sol,rjdmean_sol2] = rjd_mean(A,nb_trials);
        tmean = tmean + toc(start);

        % Best residual
        start = tic;
        rjdbestres_sol = rjd_best_res(A,normsA,nb_trials);
        tbest = tbest + toc(start);

        % Subspace RQ1
        start = tic;
        rjdsubspace_sol = rjd_subspace1(A,nb_trials);
        tsub1 = tsub1 + toc(start);
%{
        % Subspace RQ2
        start = tic;
        rjdsubspace_sol2 = rjd_subspace2(A,nb_trials);
        tsub2 = tsub2 + toc(start);
%}
        [err_mean1(k,r),~]   = residuals(problem,rjdmean_sol,@evalmon);
        [err_mean2(k,r),~]   = residuals(problem,rjdmean_sol2,@evalmon);
        [err_bestres(k,r),~] = residuals(problem,rjdbestres_sol,@evalmon);
        [err_sub1(k,r),~]    = residuals(problem,rjdsubspace_sol,@evalmon);
        %[err_sub2(k,r),~]    = residuals(problem,rjdsubspace_sol2,@evalmon);

    end

    time_mean(k)    = tmean/repeats;
    time_bestres(k) = tbest/repeats;
    time_sub1(k)    = tsub1/repeats;
    time_sub2(k)    = tsub2/repeats;

end

fprintf("\n ")
%% Export results

output_file = ['results_', example, '.tex'];
export_results(output_file, example, trial_values, repeats, ...
    err_rjd, time_rjd, ...
    err_mean1, time_mean, ...
    err_mean2, ...
    err_bestres, time_bestres, ...
    err_sub1, time_sub1)%, ...
    %err_sub2, time_sub2)

%% Print results
print_stats("Ref", err_rjd, time_rjd, repeats)
for idx = 1:length(trial_values)
    k = trial_values(idx);
    fprintf("\n   k=%d \n", k)
    print_stats('Best residual', err_bestres(idx,:), time_bestres(idx), repeats);
    print_stats('Mean RQ1', err_mean1(idx,:), time_mean(idx), repeats);
    print_stats('Mean RQ2', err_mean2(idx,:), time_mean(idx), repeats);
    print_stats('Subspace RQ1', err_sub1(idx,:), time_sub1(idx), repeats);
    print_stats('Subspace RQ2', err_sub2(idx,:), time_sub2(idx), repeats);
end

% auxiliary if want to print the stats in matlab
function print_stats(name, errors, time, rep)
    fprintf("%-15s  avg=%.2e  median=%.2e  worst=%.2e  avg_time=%.3fs\n", ...
        name, mean(errors), median(errors), max(errors), time/rep);
end
