function export_results(filename, ex_name, NB_trials, repeats, ...
    err_rjd, time_rjd, ...
    err_mean1, time_mean, ...
    err_mean2, ...
    err_bestres, time_bestres, ...
    err_sub1, time_sub1, ...
    err_sub2, time_sub2)

% Export the results of the rootfinding main to a latex table, named
% 'filename'. 

fid = fopen(filename, 'w');

% --- Begin table ---
warning off % remove warning from special chracters not existing (\c)
fprintf(fid, '\\begin{table}[!hbt!]\n \\centering \n');
fprintf(fid, ['\\caption{Rootfinding results ', ex_name, '(', num2str(repeats), ' repeats) }\n']);
fprintf(fid, '\\begin{tabular}{llcccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Method & $k$ & Avg error & Median error & Max error & Avg time (s) \\\\\n');
fprintf(fid, '\\hline\n');
% RJEA
fprintf(fid,  'Reference (RJEA) & - & %.2e & %.2e & %.2e & %.2f \\\\\n', ...
    mean(err_rjd), median(err_rjd), max(err_rjd), time_rjd);
fprintf(fid, '\\hline\n');

% My methods
print_method(fid, 'Best residual', NB_trials, err_bestres, time_bestres);
print_method(fid, 'Mean RQ1', NB_trials, err_mean1, time_mean);
print_method(fid, 'Mean RQ2', NB_trials, err_mean2, time_mean);
print_method(fid, 'Subspace RQ1', NB_trials, err_sub1, time_sub1);
if nargin >=14
    print_method(fid, 'Subspace RQ2', NB_trials, err_sub2, time_sub2);
end

% --- End table ---
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid); 

end

%% Helper to print one method with different number of trials (k)
function print_method(filename, name, trial_values, errors, times)

K = length(trial_values);
for j = 1:K

    if j == 1
        fprintf(filename, '\\multirow{%d}{*}{%s}', K, name);
    end

    fprintf(filename, ' & %d & %.2e & %.2e & %.2e & %.4f \\\\\n', ...
        trial_values(j), ...
        mean(errors(j,:)), ...
        median(errors(j,:)), ...
        max(errors(j,:)), ...
        times(j));
end
fprintf(filename, '\\hline\n');
end