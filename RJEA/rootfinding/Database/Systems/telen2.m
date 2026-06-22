function [system] = telen2()
    eqs = cell(2,1);
    eqs{1} = [1 1 0; 1/3 0 2; -1 2 0];
    eqs{2} = [-1/3 1 0; 1/3 2 0];
    
    system = systemstruct(eqs);
end