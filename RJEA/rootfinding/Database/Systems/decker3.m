function [system] = decker3() 
    eqs = cell(2,1); 
    eqs{1} = [1 1 0; 1 0 2]; 
    eqs{2} = [1.5 1 1; 1 0 2; 1 0 3]; 

    system = systemstruct(eqs); 
end