function [system] = wilkinson2(n)
    eqs = cell(1,1);
    syms x
    Q = ((x - 20)^2)*prod(x-(1:n));
    Q = expand(Q);
    Q = sym2poly(Q);
    powers = n:-1:0;
    eqs{1} = [Q' powers'];

    system = systemstruct(eqs);
end