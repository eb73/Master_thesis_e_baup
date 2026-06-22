function [system] = telen4()
    eqs = cell(2,1);
    eqs{1} = [1 0 0; 1 1 0; 1 0 1; 1 1 1; 1 2 1; 1 3 1];
    eqs{2} = [1 0 0; 1 0 1; 1 1 1; 1 2 1];

    system = systemstruct(eqs);
end