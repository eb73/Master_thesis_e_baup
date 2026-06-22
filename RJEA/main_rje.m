%% main_rje.m
% Main file used for the numerical experiments.
% Relies on a struct 'user' or default variables, to specify setting.
% user may contain: example (which matrices to use), Nsamples (number of
% runs), nb_trials (number of initial trials used = k), rq_type
% ('one-side', 'two-sided', 'both'), method

addpath 'RJEA_He_Plestenjak'  % for reference method & other auxiliaries

%% Setting

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

% verify that parameters are supported by chosen method
supports_2_sided = {'ref', 'smaller_subspace', 'best_res', 'mean'}; 
needs_two_sided = ismember(rq_type, {'two-sided', 'both'}); % flag rq_type that needs RQ2

if needs_two_sided && ~ismember(method, supports_2_sided)
    warning('%s method does not support two-sided RQ. Falling back to one-sided.', method)
    rq_type = 'one-sided';
end

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
    str_met_clean = strcat(num2str(S.nb_trials), 'D subspace');
end

if ismember(method, {'best_res', 'mean'})
    str_met_clean = strcat(str_met_clean, ' (k=', num2str(S.nb_trials), ')');
end

str_output = strcat('ex', str_ex, '_', str_met);

fprintf('Running example %3.1f with method %s and Nsamples=%d',S.example,str_met_clean,Nsamples)
if strcmp('ref', method)~=1
    str_output = strcat(str_output, num2str(S.nb_trials));
end
fprintf(', RQ type: %s', rq_type)
fprintf('\n')

materr2SA = [];
materr1SA = [];
materr2X = [];
materr1X = [];
cndhist = [];
cndhist2 = [];
res1A = [];
res2A = [];
bound_simple1 = [];
bound_simple2_eps = [];
bound_simple2 = [];
bound_schur1 = [];
bound_schur2 = [];
nX1A = [];
nX2A = [];
nY1A = [];
nY2A = [];

%% Run experiments 
for ind = 1:length(noise_level)
    fprintf('computing for noise %5.2e \n',noise(ind))
    A = S.Aset;
    for j = 1:S.d 
        A{j} = A{j} + noise(ind)/sqrt(S.d)*S.PA{j};
    end
    use_noise = max(noise(ind),S.min_noise);
    %[rq1, res, err, cnd, d2, normsX] = newRayleighSmallerSubspace(A, Nsamples, S.exact_eig, S.index, 1, S.noise_for_Amu, S.Aset, S.nb_trials);

    [rq1, rq2, res1, err1, res2, err2, cnd, d2, normsX] = newRayleighOpt(A, S.Aset, Nsamples, 1, S);
    [~,ord2] = sort(rownorms(err2(:,1)));
    [~,ord1] = sort(rownorms(err1(:,1)));
    materr1SA = [materr1SA sort(rownorms(err1(:,1)))];
    materr1X = [materr1X sort(rownorms(err1(:,2)))];
    materr2SA = [materr2SA sort(rownorms(err2(:,1)))];
    materr2X = [materr2X sort(rownorms(err2(:,2)))];
    cndhist = [cndhist sort(cnd(:,1))];
    if size(cnd, 2) == 2 % case smaller_subspace+both
        cndhist2 = [cndhist2 sort(cnd(:,2))];
    end
    boundRQ1 = normsX(:,1).^2.*(S.delta + use_noise*(1 + sqrt(2)*normsX(:,2).*normsX(:,4).*d2(:,1)));
    bound_simple1 = [bound_simple1 boundRQ1(ord1)];
    boundRQ2_eps = use_noise*normsX(:,1).*normsX(:,3);
    boundRQ2 = boundRQ2_eps.*(1 + use_noise*normsX(:,2).*normsX(:,4).*(sqrt(2)*d2(:,3)));
    bound_simple2 = [bound_simple2 boundRQ2(ord2)];
    res1A = [res1A sort(res1(:,1))];
    res2A = [res2A sort(res2(:,1))];
    nX1A = [nX1A sort(normsX(:,1))]; 
    nX2A = [nX2A sort(normsX(:,2))];
    nY1A = [nY1A sort(normsX(:,3))];
    nY2A = [nY2A sort(normsX(:,4))];
    
end

%% Plots

legend_text = cell(1,length(noise_level));
for j = 1:length(noise_level)
    legend_text{j} = sprintf('10^{%d}',noise_level(j));
end
if isinf(noise_level(1))
    legend_text{1} = 'no noise';
end

rq_out = '';
if plot_onesided
    rq_out = '_rq1';
elseif plot_twosided
    rq_out = '_rq2';
end

% Figure 1: errors

% Error distribution of 1-S
if find(plot_figures==1) && plot_onesided
    f1 = figure;
    Ax1(1) = axes(f1);
    h = semilogy(materr1X, '--');
    lgd1 = legend(h, legend_text,'Location','southeast');
    title(sprintf('Error distribution of RQ1 for \\lambda= %s, %s',string_eig, str_met_clean)) 
    ylim(S.ylim_err)
    % export if needed
    if download_output
        exportgraphics(f1,strcat('figures\', str_output, '_err', rq_out, '.png'))
    end
end

% Error distribution of 2-S
if find(plot_figures==1) && plot_twosided
    f1 = figure;
    Ax1(1) = axes(f1);
    h = semilogy(materr2X);
    lgd1 = legend(h, legend_text,'Location','southeast');
    title(sprintf('Error distribution of RQ2 for \\lambda= %s, %s',string_eig, str_met_clean)) 
    ylim(S.ylim_err)
    % export if needed
    if download_output
        exportgraphics(f1,strcat('figures\', str_output, '_err', rq_out, '.png'))
    end
end

% Error distribution of 1-S & 2-S in same graph
if find(plot_figures==1) && plot_both
    f1 = figure;
    Ax1(1) = axes(f1);
    semilogy(materr2X(:,1))
    hold on
    for ind = 2:length(noise_level)
        semilogy(materr2X(:,ind))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        semilogy(materr1X(:,ind),'--')
    end
    lgd1 = legend(legend_text,'Location','southeast');
    title(sprintf('Error distribution for \\lambda= %s, %s',string_eig, str_met_clean))
    ylim(S.ylim_err)
    Ax1(2) = copyobj(Ax1(1),gcf);
    delete(get(Ax1(2),'children'))
    
    % below is to plot legend of RQ1 vs RQ2
    hold on
    H1 = plot(nan, nan, '--', 'Color', [0 0 0], 'Parent', Ax1(2), 'Visible', 'on');
    H2 = plot(nan, nan, '-', 'Color', [0 0 0], 'Parent', Ax1(2), 'Visible', 'on');
    hold off
    set(Ax1(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd2 = legend([H1 H2], '1-sided RQ', '2-sided RQ', 'Location', 'south');
    set(lgd2,'color','white')
    set(lgd2, 'TextColor', 'black');
    
    % export if needed
    if download_output
        exportgraphics(f1,strcat('figures\', str_output, '_err.png'))
    end
end

% Figure 2: residuals

if find(plot_figures==2) && plot_onesided
    f2 = figure;
    plot(log10(res1A(:,1)), '--')
    hold on
    for ind = 2:length(noise_level)
        plot(log10(res1A(:,ind)), '--')
    end
    legend(legend_text, 'Location', 'southeast')
    title(sprintf('Residual distribution RQ1, %s', str_met_clean))
    ylim(S.ylim_res)
    % export
    if download_output
        exportgraphics(f2,strcat('figures\', str_output, '_resnorm', rq_out, '.png'))
    end
end

if find(plot_figures==2) && plot_twosided
    f2 = figure;
    plot(log10(res2A(:,1)))
    hold on
    for ind = 2:length(noise_level)
        plot(log10(res2A(:,ind)))
    end
    legend(legend_text, 'Location', 'southeast')
    title(sprintf('Residual distribution RQ2, %s', str_met_clean))
    ylim(S.ylim_res)
    % export
    if download_output
        exportgraphics(f2,strcat('figures\', str_output, '_resnorm', rq_out, '.png'))
    end
end

if find(plot_figures==2) && plot_both
    f2 = figure;
    Ax2 = axes(f2);
    plot(log10(res2A(:, 1)))
    hold on
    for ind = 2:length(noise_level)
        plot(log10(res2A(:,ind)))
    end
    set(gca,'ColorOrderIndex',1)
    for ind = 1:length(noise_level)
        plot(log10(res1A(:, ind)), '--')
    end
    lg2 = legend(legend_text, 'Location', 'southeast');
    title(sprintf('Residual distribution (log10 scale), %s', str_met_clean))
    ylim(S.ylim_res)
    Ax2(2) = copyobj(Ax2(1),gcf);
    delete(get(Ax2(2),'children'))
    
    % below is to plot legend of RQ1 vs RQ2
    hold on
    H1 = plot(nan, nan, '--', 'Color', [0 0 0], 'Parent', Ax2(2), 'Visible', 'on');
    H2 = plot(nan, nan, '-', 'Color', [0 0 0], 'Parent', Ax2(2), 'Visible', 'on');
    hold off
    set(Ax2(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')
    lgd2 = legend([H1 H2], '1-sided RQ', '2-sided RQ', 'Location', 'south');
    set(lgd2,'color','white')
    set(lgd2, 'TextColor', 'black');
    % export
    if download_output
        exportgraphics(f2,strcat('figures\', str_output, '_resnorm.png'))
    end
end

%% Tables

% Build data 
Med = round(Nsamples/2);
switch rq_type
    case 'one-sided'
        row_median = [noise; materr1X(Med,:); res1A(Med,:); cndhist(Med,:)]';
        row_worst  = [noise; materr1X(end,:); res1A(end,:); cndhist(end,:)]';
        hdr_short  = 'eps      | 1S RQ err | 1S RQ res |  cond(X)  |';
        colnames   = {'noise', 'error', 'residual', 'condX'};
        col_idx    = 1:4;

    case 'two-sided'
        row_median = [noise; materr2X(Med,:); res2A(Med,:); cndhist(Med,:)]';
        row_worst  = [noise; materr2X(end,:); res2A(end,:); cndhist(end,:)]';
        hdr_short  = 'eps      | 2S RQ err | 2S RQ res |  cond(X)  |';
        colnames   = {'noise', 'error', 'residual', 'condX'};
        col_idx    = 1:4;

    case 'both'
        if strcmp(method,'smaller_subspace')
            row_median = [noise; materr1X(Med,:); res1A(Med,:); cndhist(Med,:); materr2X(Med,:); res2A(Med,:); cndhist2(Med,:)]';
            row_worst  = [noise; materr1X(end,:); res1A(end,:); cndhist(Med,:); materr2X(end,:); res2A(end,:); cndhist2(end,:)]';
            hdr_short  = 'eps      | 1S RQ err | 1S RQ res |1S RQ cond(X)| 2S RQ err | 2S RQ res |2S RQ cond(X)|';
            colnames   = {'noise', '1S error', '1S residual', '1S condX', '2S error', '2S residual', '2S condX'};
            col_idx    = 1:7;
        else

            row_median = [noise; materr1X(Med,:); res1A(Med,:); materr2X(Med,:); res2A(Med,:); cndhist(Med,:)]';
            row_worst  = [noise; materr1X(end,:); res1A(end,:); materr2X(end,:); res2A(end,:); cndhist(end,:)]';
            hdr_short  = 'eps      | 1S RQ err | 1S RQ res | 2S RQ err | 2S RQ res |  cond(X)  |';
            colnames   = {'noise', '1S error', '1S residual', '2S error', '2S residual', 'condX'};
            col_idx    = 1:6;
        end
end

divider = [repmat('-', 1, length(hdr_short)-1), '|'];

% Print
fprintf('\nMedian error and condition numbers\n')
fprintf('%s\n', hdr_short);
fprintf('%s\n', divider);
for k = 1:size(row_median, 1)
    fmt = repmat('  %7.1e  |', 1, length(col_idx));
    fprintf(['%5.1e  |', fmt], row_median(k, col_idx));
    fprintf('\n')
end

fprintf('\nWorst error and condition numbers\n')
fprintf('%s\n', hdr_short);
fprintf('%s\n', divider);
for k = 1:size(row_worst, 1)
    fmt = repmat('  %7.1e  |', 1, length(col_idx));
    fprintf(['%5.1e  |', fmt], row_worst(k, col_idx));
    fprintf('\n')
end

% Export
if download_output
    if plot_both
        latex_table(strcat('tables/', str_output, '_median', '.tex'), strcat('Median results for example ', ' ', str_ex, ', method: ', str_met_clean, ' (', num2str(Nsamples), ' trials)'), row_median(:, col_idx), colnames)
        latex_table(strcat('tables/', str_output, '_worst', '.tex'), strcat('Worst results for example ', ' ', str_ex, ', method: ', str_met_clean, ' (', num2str(Nsamples), ' trials)'), row_worst(:, col_idx), colnames)
    else
        if plot_onesided
            txt_result = [' results RQ1 for example', ' ', str_ex, ', method: ', str_met_clean, ' (', num2str(Nsamples), ' trials)'];
            latex_table(strcat('tables/', str_output, '_median', rq_out, '.tex'), strcat('Median', txt_result), row_median(:, col_idx), colnames)
            latex_table(strcat('tables/', str_output, '_worst', rq_out, '.tex'), strcat('Worst', txt_result), row_worst(:, col_idx), colnames)
        else
            txt_result = [' results RQ2 for example', ' ', str_ex, ', method: ', str_met_clean, ' (', num2str(Nsamples), ' trials)'];
            latex_table(strcat('tables/', str_output, '_median', rq_out, '.tex'), strcat('Median', txt_result), row_median(:, col_idx), colnames)
            latex_table(strcat('tables/', str_output, '_worst', rq_out, '.tex'), strcat('Worst', txt_result), row_worst(:, col_idx), colnames)
        end
        colnames_combined = {'noise', 'median error', 'median residual', 'median cond(X)', 'worst error', 'worst residual', 'worst cond(X)'};
        txt_result = ['Results for example ', str_ex, ', method: ', str_met_clean, ' (', num2str(Nsamples), ' trials)'];
        latex_table(strcat('tables/', str_output, rq_out, '.tex'), txt_result, [row_median(:, col_idx), row_worst(:, 2:4)], colnames_combined)
    end
end
