function latex_table(filename, table_caption, data, colnames)
% write the data as a latex table in the file filename
% filename     : (string) with name/path of the tex file to write the table in
% tablecaption : (string) caption to put for the table
% data         : (n x m) numeric matrix
% colnames     : 1 x m cell array of column names

[n, m] = size(data);

% verify that colnames is the same size as #cols of data
if m~=length(colnames)
    error("Data and names of columns don't have the same size")
end

fid = fopen(filename, 'w'); % open file if exists, creates it otherwise

% --- Begin table ---
warning off % remove warning from special chracters not existing (\c)
fprintf(fid, '\\begin{table}[!hbt!]\n \\centering \n');
fprintf(fid, strcat('\\caption{', table_caption, '}\n'));
fprintf(fid, '\\begin{tabular}{');

fprintf(fid, repmat('c',1,m));
fprintf(fid, '}\n\\hline\n');
warning on
% --- Header ---

for j = 1:m
    if j == m
        fprintf(fid, '%s \\\\\n', colnames{j});
    else
        fprintf(fid, '%s & ', colnames{j});
    end
end
fprintf(fid, '\\hline\n');

% --- Data ---
for i = 1:n
    
    for j = 1:m
        if j == 1
            val = latex_format(data(i,j), true); % header is noise level
        else
            val = latex_format(data(i,j));
        end
        if j == m
            fprintf(fid, '%s \\\\\n', val); % last line: add \\ for newline
        else
            fprintf(fid, '%s & ', val); % other lines: val and tab symbol to make table
        end
    end
end

% --- End table ---
fprintf(fid, '\\hline\n\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid); % close the file to ensure good behavior at the end
end

%% Helper function for latex formatting
function str = latex_format(x, power10)
    if nargin<2, power10 = false; end

    if x == 0
        str = '0';
        return;
    end
    
    if power10
        e = round(log10(abs(x)));
        str = sprintf('$10^{%d}$', e);
    else
        e = floor(log10(abs(x)));
        m = x / 10^e;
        % avoid cases where mantissa to 1 decimal gives 10.0
        if round(m, 1) >= 10
            e = e + 1;
            m = x / 10^e;
        end
        % format mantissa
        if e==0
            % don't write 10^0
            str = sprintf('$%.1f$', m);
        else
            str = sprintf('$%.1f \\cdot 10^{%d}$', m, e);
        end
    end
end
