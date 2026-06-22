function [system] = bronstein(a,b)
    eqs = cell(3,1);
    eqs{1} = [2 0 4; 2 2 2; (b^2-3*a^2) 0 2; -2*b 1 2; -2*b 0 3; 2*a^2*b 1 0; 2*a^2*b 0 1; -a^2 2 0; a^2*(a^2-b^2) 0 0];
    eqs{2} = [8 0 3; 4 2 1; -6*b 0 2; -4*b 1 1; 2*(b^2-3*a^2) 0 1; 2*a^2*b 0 0];
    eqs{3} = [4 1 2; -2*b 0 2; -2*a^2 1 0; 2*a^2*b 0 0];

    system = systemstruct(eqs);
end