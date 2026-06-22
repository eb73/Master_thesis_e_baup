function [system] = batselier8() 
    eqs = cell(3,1); 
    eqs{1} = [1 2 0 0; 1 1 0 1; -2 0 1 0; 5 0 0 0]; 
    eqs{2} = [2 3 1 0; 7 0 1 2; -4 1 1 1; 3 1 0 0; -2 0 0 0]; 
    eqs{3} = [1 0 4 0; 2 0 1 1; 5 2 0 0; -5 0 0 0]; 

    system = systemstruct(eqs); 
end