function [system] = batselier11() 
    eqs = cell(2,1); 
    eqs{1} = [1 0 2; -6 0 1; 9 0 0]; 
    eqs{2} = [1 2 0; -2 1 1; 2 1 0; 1 0 2; -2 0 1; 1 0 0]; 

    system = systemstruct(eqs); 
end