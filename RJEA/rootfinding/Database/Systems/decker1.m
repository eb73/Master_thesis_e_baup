function [system] = decker1() 
    eqs = cell(2,1); 
    eqs{1} = [1 3 0; 1 1 1]; 
    eqs{2} = [1 0 2; 1 0 1]; 

    system = systemstruct(eqs); 
end