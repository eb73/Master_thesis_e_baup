function [system] = large2(n,a)
    eqs = cell(1,1);
    eqs{1} = [a n; 1e300 1; 1 11; 1e-300 0];

    system = systemstruct(eqs);
end