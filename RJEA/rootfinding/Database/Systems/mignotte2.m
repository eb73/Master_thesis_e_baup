function [system] = mignotte2(n,m,a,k)
    eqs = cell(1,1);
    C = zeros(1,m+1);
    D = zeros(1,k+1);
    powers = 0:m+k;
    for l = 0:m+k
        C(l+1) = nchoosek(m+k,l)*a^l;
    end
    powers2 = n-k:n;
    for l = 0:k
        D(l+1) = nchoosek(k,l)*a^l;
    end
    eqs{1} = [C' powers'; D' powers2'];

    system = systemstruct(eqs);
end