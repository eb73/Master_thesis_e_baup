function [system] = mickey() 
    eqs = cell(2,1); 
    eqs{1} = [1 2 0; 4 0 2; -4 0 0]; 
    eqs{2} = [-1 1 0; 2 0 2]; 

    system = systemstruct(eqs); 
end