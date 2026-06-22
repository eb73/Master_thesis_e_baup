%% analysis_n_ex1.m
% File for analysis of methods depending on matrices size (n)
% The struct user is overridden

addpath 'RJEA_He_Plestenjak'   % adjust path if needed

%% Setting parameters
n_values   = [5,10];   % matrix sizes to check (multiples of 5)
Nsamples   = 100;       
noise_idx  = 2;                       % which noise level to report (index into noise_level array)

methods = {'ref', 'best_res', 'mean', 'smaller_subspace'};

% For each n and method: median & worst error and residual
med_err1  = nan(length(n_values), length(methods));
wst_err1  = nan(length(n_values), length(methods));
med_res1  = nan(length(n_values), length(methods));
wst_res1  = nan(length(n_values), length(methods));

med_err2  = nan(length(n_values), length(methods));
wst_err2  = nan(length(n_values), length(methods));
med_res2  = nan(length(n_values), length(methods));
wst_res2  = nan(length(n_values), length(methods));

%% Run experiment
for im = 1:length(methods)
    opt = methods{im};
    use_two = ~strcmp(opt, 'smaller_subspace');  % subspace: one-sided only

    for in = 1:length(n_values)
        n_cur = n_values(in);
        fprintf('[%s]  n = %d ...\n', opt, n_cur);

        % Build example 1.1 first to setup all variables
        user.example = 1.1; % base example 1 (n=10 default)
        user.method = opt;
        user.Nsamples = Nsamples;
        user.download_output = false;
        user.rq_type = 'both'; 
        % for subspace method, don't show rq2 since not sound
        if ~use_two
            user.rq_type = 'one-sided';
        end
        S = config(user);

        % override matrices for different n
        rng(1)
        da_n = repelem((1:n_cur/5), 5)';
        db_n = repmat((1:5), 1, n_cur/5)';
        DA_n = diag(da_n);
        DB_n = diag(db_n);
        kappa = 1e2;
        class_t = 'double';
        X1_n = condmatX(n_cur, kappa, class_t);
        Y1_n = inv(X1_n);
        A1_n = X1_n * DA_n * Y1_n;
        B1_n = X1_n * DB_n * Y1_n;

        S.A1 = A1_n;
        S.B1 = B1_n;
        S.Aset = {A1_n, B1_n};
        S.exact_eig = [da_n db_n];
        S.index = 1;
        S.d = 2;
        S.delta = 0;

        % Rebuild noise matrices for new size
        realfamily = isreal(A1_n) && isreal(B1_n);
        PA_n = cell(2,1);
        for j = 1:2
            if realfamily
                PA_n{j} = randn(n_cur, class_t);
            else
                PA_n{j} = randn(n_cur, class_t) + 1i*randn(n_cur, class_t);
            end
            PA_n{j} = PA_n{j} / norm(PA_n{j}, 'fro');
        end
        S.PA       = PA_n;
        S.min_noise = norm([norm(A1_n) norm(B1_n)]) * eps(class_t) / 2;
        gapA_n = Inf;
        for j = 2:n_cur
            gapA_n = min(norm(da_n(j,:) - da_n(1,:)), gapA_n); 
        end
        S.gapA = gapA_n;

        % pick noise level
        noise_level_arr = S.noise_level;
        noise_arr       = S.noise;
        ni = min(noise_idx, numel(noise_arr));

        A_run = S.Aset;
        for j = 1:S.d
            A_run{j} = A_run{j} + noise_arr(ni)/sqrt(S.d) * S.PA{j};
        end
        use_noise_val = max(noise_arr(ni), S.min_noise);

        % Run method
        [rq1, rq2, res1, err1, res2, err2, cnd, d2, normsX] = newRayleighOpt(A_run, S.Aset, Nsamples, 1, S);

        Med = round(Nsamples / 2);

        % one-sided
        e1 = sort(rownorms(err1(:,2)));  
        r1 = sort(res1(:,1));             
        med_err1(in, im) = e1(Med);
        wst_err1(in, im) = e1(end);
        med_res1(in, im) = r1(Med);
        wst_res1(in, im) = r1(end);

        % two-sided (only when available)
        if use_two
            e2 = sort(rownorms(err2(:,2)));
            r2 = sort(res2(:,1));
            med_err2(in, im) = e2(Med);
            wst_err2(in, im) = e2(end);
            med_res2(in, im) = r2(Med);
            wst_res2(in, im) = r2(end);
        end
    end
end


%% Plots
% Setting
colors  = lines(length(methods));
ms = 7;
nb_trials = S.nb_trials;
noise_level = S.noise_level(noise_idx);
method_labels = {'ref', sprintf('best\\_res (k=%d)',nb_trials), ...
                         sprintf('mean (k=%d)',nb_trials), ...
                         sprintf('subspace (k=%d)',nb_trials)};
ylabels     = {'error', 'error', 'residual', 'residual'};

%% Figure 1: RQ1, absolute errors & residuals
f1 = figure('Name','Errors & residuals vs n (RQ1)','Position',[50 50 1100 820]);
tl1  = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');
title(tl1, sprintf('Errors & residuals vs n (RQ1), noise = 10^{%d}', ...
    noise_level), 'FontSize', 12)

panel_data  = {med_err1, wst_err1, med_res1, wst_res1};
panel_title = {'Median error  (RQ1)', 'Worst error  (RQ1)', ...
               'Median residual  (RQ1)', 'Worst residual  (RQ1)'};

for ip = 1:4
    nexttile(tl1)
    hold on; box on; grid on
    for im = 1:length(methods)
        semilogy(n_values, panel_data{ip}(:,im), ...
            '-', 'Color', colors(im,:), 'Marker', 'x', ...
            'LineWidth', 1, 'MarkerSize', ms, 'DisplayName', method_labels{im})
    end
    xlabel('n'); ylabel(ylabels{ip})
    title(panel_title{ip})
    legend('Location','best','Interpreter','tex','FontSize',8)
    set(gca,'YScale','log','XTick',n_values)
end
if S.download_output
    exportgraphics(f1,strcat('ex1.1_n_analysis_rq1_k', num2str(nb_trials), '_', num2str(noise_level), '.png'))
end

%% Figure 2: RQ2, absolute errors & residuals 
f2 = figure('Name','Errors & residuals vs n (RQ2)','Position',[50 50 1100 820]);
tl2  = tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');
title(tl2, sprintf('Errors & residuals vs n (RQ2), noise = 10^{%d}', ...
    noise_level), 'FontSize', 12)

panel_data  = {med_err2, wst_err2, med_res2, wst_res2};
panel_title = {'Median error  (RQ2)', 'Worst error  (RQ2)', ...
               'Median residual  (RQ2)', 'Worst residual  (RQ2)'};

for ip = 1:4
    nexttile(tl2)
    hold on; box on; grid on
    for im = 1:length(methods)-1
        semilogy(n_values, panel_data{ip}(:,im), ...
            '-', 'Color', colors(im,:), 'Marker', 'x', ...
            'LineWidth', 1, 'MarkerSize', ms, 'DisplayName', method_labels{im})
    end
    xlabel('n'); ylabel(ylabels{ip})
    title(panel_title{ip})
    legend('Location','best','Interpreter','tex','FontSize',8)
    set(gca,'YScale','log','XTick',n_values)
end

if S.download_output
    exportgraphics(f2,strcat('ex1.1_n_analysis_rq2_k', num2str(nb_trials), '_', num2str(noise_level), '.png'))
end

fprintf('\nAnalysis done.\n')