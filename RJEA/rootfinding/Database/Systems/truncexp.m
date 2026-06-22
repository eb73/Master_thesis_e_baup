function [system] = truncexp(m)
    eqs = cell(1,1);
    M = zeros(m+1,2);
    for k = 0:m
        M(k+1,1) = 1/factorial(k);
        M(k+1,2) = k;
    end
    eqs{1} = M;

    system = systemstruct(eqs);
end