function [system] = kameny1(c)
    eqs = cell(1,1);
    eqs{1} = [9*c^4 0; -6*c^2 1; 1i*c 7; 1 2];

    system = systemstruct(eqs);
end