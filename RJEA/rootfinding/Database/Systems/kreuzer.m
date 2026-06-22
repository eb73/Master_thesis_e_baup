function [system] = kreuzer() 
    eqs = cell(2,1); 
    eqs{1} = [0.25 2 0; 1 0 2; -1 0 0]; 
    eqs{2} = [1 2 0; 0.25 0 2; -1 0 0]; 

    system = systemstruct(eqs); 
end