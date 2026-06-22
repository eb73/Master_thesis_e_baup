function [system] = phcpack1() 
    eqs = cell(2,1); 
    eqs{1} = [1 2 0]; 
    eqs{2} = [1 0 2]; 

    system = systemstruct(eqs); 
end