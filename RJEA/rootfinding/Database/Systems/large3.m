function [system] = large3(n,a)
    eqs = cell(1,1);
    eqs{1} = [a*1e-200 n; 1e100 n-1; 1e200 0];

    system = systemstruct(eqs);
end