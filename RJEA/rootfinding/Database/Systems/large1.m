function [system] = large1(n,a)
    eqs = cell(1,1);
    eqs{1} = [a n; 1e300 14; 1 5; 1 0];

    system = systemstruct(eqs);
end