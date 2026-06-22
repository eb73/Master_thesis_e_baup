function [system] = nthroot(n,a)
    eqs = cell(1,1);
    eqs{1} = [1 n; -a 0];

    system = systemstruct(eqs);
end