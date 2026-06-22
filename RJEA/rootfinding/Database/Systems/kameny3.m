function [system] = kameny3(c)
    eqs = cell(1,1);
    eqs{1} = [c^4 4; c^2 9; -6*c^2 2; 9 0];

    system = systemstruct(eqs);
end