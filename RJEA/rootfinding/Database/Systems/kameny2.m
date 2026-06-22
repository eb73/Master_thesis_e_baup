function [system] = kameny2(c)
    eqs = cell(1,1);
    eqs{1} = [c^4 4; 1i*c^2 9; -6*c^2 2; 9 0];

    system = systemstruct(eqs);s
end