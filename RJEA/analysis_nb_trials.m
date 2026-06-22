%% analysis_nb_trials.m 
% File for analysis of methods depending on nb_trials (k)
% Relies on same struct that main_rje.m, field nb_trials is overridden

addpath 'RJEA_He_Plestenjak' % for reference method & auxiliaries

%% Setting

NB_trials = 1:6; % change if wanna investigate different number of trials (k)

if exist('user', 'var')
    S = config(user);
else
    S = config();
end

% these are parameters that are used a few times, so avoid S. each time for readability
plot_figures = S.plot_figures;
string_eig = S.string_eig;
noise_level = S.noise_level;
noise = S.noise; 
method = S.method; 
Nsamples = S.Nsamples;
download_output = S.download_output;
rq_type = S.rq_type;

% some variables to help with tables and figures
plot_onesided = strcmp(rq_type, 'one-sided');
plot_twosided = strcmp(rq_type,'two-sided');
plot_both = strcmp(rq_type, 'both');
str_ex = num2str(S.example);
str_met = method;
str_met_clean = str_met;
if strcmp(method, 'best_res')
    str_met_clean = 'best residual';
elseif strcmp(method, 'smaller_subspace')
    str_met_clean = 'smaller subspace';
end

str_output = strcat('ex', str_ex, '_', str_met);

fprintf('Running analysis of method vs number of trials (k) for example %3.1f with %s and Nsamples=%d',S.example,str_met_clean,Nsamples)
fprintf(', RQ type: %s', rq_type)
fprintf('\n')

rq_out = '';
if plot_onesided
    rq_out = '_rq1';
elseif plot_twosided
    rq_out = '_rq2';
end

materr2SA = [];
materr1SA = [];
materr2X = [];
materr1X = [];
cndhist = [];
res1A = [];

median_res = zeros(length(noise_level), length(NB_trials));
worst_res = zeros(length(noise_level), length(NB_trials));
median_err = zeros(length(noise_level), length(NB_trials));
worst_err = zeros(length(noise_level), length(NB_trials));
median_idx = round(Nsamples/2);

nb_repet = 10; % number of repetitions

median_res_all = zeros(length(noise_level), length(NB_trials), nb_repet);
median_err_all = zeros(length(noise_level), length(NB_trials), nb_repet);
worst_res_all  = zeros(length(noise_level), length(NB_trials), nb_repet);
worst_err_all  = zeros(length(noise_level), length(NB_trials), nb_repet);

%% Run method
for rep = 1:nb_repet
    fprintf('Repetition %d/%d\n', rep, nb_repet);
    for idx = 1:length(noise_level)

        fprintf('computing for noise %5.2e \n',noise(idx))
        A = S.Aset;
        for j = 1:S.d
            A{j} = A{j} + noise(idx)/sqrt(S.d)*S.PA{j};
        end
        for idx2 = 1:length(NB_trials)
            nb_trials = NB_trials(idx2);
            if nb_trials == 1
                % 1 trial corresponds to reference method
                fprintf('running for 1 trial (reference)\n')
                [~, rq1, res, err, ~, ~, ~, ~] = testRayleigh(A,Nsamples,S.exact_eig,S.index,1,S.noise_for_Amu,0,S.Aset,0);
                materr1SA = sort(rownorms(err(:,2)));
                materr1X = sort(rownorms(err(:,4)));
                res1A = sort(res(:,2));
            else
                S.nb_trials = nb_trials;
                fprintf('running for %d trials \n', S.nb_trials)
                
                [rq1, rq2, res1, err1, res2, err2, cnd, d2, normsX] = newRayleighOpt(A, S.Aset, Nsamples, 1, S);
                materr1SA = sort(rownorms(err1(:,1)));
                materr1X = sort(rownorms(err1(:,2)));
                res1A = sort(res1(:,1));
            end

            median_res(idx, idx2) = res1A(median_idx);
            median_err(idx, idx2) = materr1X(median_idx);
            worst_res(idx, idx2) = res1A(end);
            worst_err(idx, idx2) = materr1X(end);
            median_res_all(idx, idx2, rep) = res1A(median_idx);
            median_err_all(idx, idx2, rep) = materr1X(median_idx);
            worst_res_all(idx, idx2, rep)  = res1A(end);
            worst_err_all(idx, idx2, rep)  = materr1X(end);
        end
    end
end


%% Plots

for idx = 1:length(noise_level)
    f = figure('Position', [100, 100, 850, 400]);

    % Residuals 
    subplot(1,2,1)
    semilogy(NB_trials, median_res(idx,:), '-o', 'DisplayName','Median')
    hold on
    semilogy(NB_trials, worst_res(idx,:), '-s', 'DisplayName','Worst')
    hold off
    title('Residuals')
    xlabel('Number of trials')
    ylabel('Residual')
    legend('Location','best')
    grid on

    % Errors
    subplot(1,2,2)
    semilogy(NB_trials, median_err(idx,:), '-o', 'DisplayName','Median')
    hold on
    semilogy(NB_trials, worst_err(idx,:), '-s', 'DisplayName','Worst')
    hold off
    title('Errors')
    xlabel('Number of trials')
    ylabel('Error')
    legend('Location','best')
    grid on

    if idx == 1
        sgtitle('Results for no noise')
    else
        sgtitle(sprintf('Results for noise = 10^{%d}', noise_level(idx)))
    end
    if download_output
        exportgraphics(f,strcat('figures\', str_output, rq_out, '_nbtrials_', num2str(noise_level(idx)),  '.png'))
    end
end


% variability of results
mean_median_res = mean(median_res_all, 3);
std_median_res  = std(median_res_all, 0, 3);

mean_median_err = mean(median_err_all, 3);
std_median_err  = std(median_err_all, 0, 3);

mean_worst_res = mean(worst_res_all, 3);
std_worst_res  = std(worst_res_all, 0, 3);

mean_worst_err = mean(worst_err_all, 3);
std_worst_err  = std(worst_err_all, 0, 3);

for idx = 1:length(noise_level)
    f1 = figure;
    hold on
    x = NB_trials;
    
    y = mean_median_res(idx,:);
    e = std_median_res(idx,:);
    fill([x fliplr(x)], [y-e fliplr(y+e)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
    semilogy(x, y, 'b-o', 'DisplayName','Median')
    
    y = mean_worst_res(idx,:);
    e = std_worst_res(idx,:);
    fill([x fliplr(x)], [y-e fliplr(y+e)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility', 'off')
    semilogy(x, y, 'r-s', 'DisplayName','Worst')
    
    legend
    if idx == 1
        title('Residuals (mean +/- std), no noise')
    else
        title(sprintf('Residuals (mean +/- std), noise = 10^{%d}', noise_level(idx)))
    end
    grid on
    hold off
    if download_output
        exportgraphics(f1,strcat('figures\', str_output, rq_out, '_variability_res_nbtrials_', num2str(noise_level(idx)),  '.png'))
    end
    f2 = figure;
    errorbar(NB_trials, mean_median_res(idx,:), std_median_res(idx,:), '-o')
    if idx == 1
        title('Errors, no noise')
    else
        title(sprintf('Errors, noise = 10^{%d}', noise_level(idx)))
    end
    if download_output
        exportgraphics(f2,strcat('figures\', str_output, rq_out, '_variability_err_nbtrials_', num2str(noise_level(idx)),  '.png'))
    end
end
