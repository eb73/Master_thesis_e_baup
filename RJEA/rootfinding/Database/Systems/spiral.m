function [system] = spiral(n,a)
    eqs = cell(1,1);
    syms x
    Q = prod(x+((a.^(1:n)-1)/(a-1)));
    Q = expand(Q);
    Q = sym2poly(Q);
    powers = n:-1:0;
    eqs{1} = [Q' powers'];

    system = systemstruct(eqs);
end