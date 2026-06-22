function [system] = ojika1() 
    eqs = cell(2,1); 
    eqs{1} = [1 2 0; 1 0 1; -3 0 0]; 
    eqs{2} = [1 1 0; 0.125 0 2; -1.5 0 0]; 

    system = systemstruct(eqs); 
end