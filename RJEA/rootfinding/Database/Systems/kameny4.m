function [system] = kameny4()
    eqs = cell(1,1);
    eqs{1} = [1 14; 2e24 11; 1e48 8; 4 7; -4e24 4; 4 0];

    system = systemstruct(eqs);
end