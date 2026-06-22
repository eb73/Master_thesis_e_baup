function [system] = telen3()
    eqs = cell(2,1);
    eqs{1} = [1 2 0; -3 1 1; 2 0 2; 1 0 0];
    eqs{2} = [1 2 0; -1 0 2; -3 0 1; 1 0 0];

    system = systemstruct(eqs);
end