function [system] = mignotte1(n,m,a)
    eqs = cell(1,1);
    C = zeros(m+1,1);
    powers = 0:m;
    for k = 0:m
        C(k+1) = nchoosek(m,k)*a^k;
    end
    eqs{1} = [C' powers'; 1 n];

    system = systemstruct(eqs);
end