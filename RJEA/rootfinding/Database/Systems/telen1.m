function [system] = telen1()
    eqs = cell(2,1);
    eqs{1} = [7 0 0; 3 1 0; -6 0 1; -4 2 0; 2 1 1; 5 0 2];
    eqs{2} = [-1 0 0; -3 1 0; 14 0 1; -2 2 0; 2 1 1; -3 0 2];

    system = systemstruct(eqs);
end