function [system] = batselier5() 
    eqs = cell(3,1); 
    eqs{1} = [1 12 0 0; 1 0 12 0; 1 0 0 12; -4 0 0 0]; 
    eqs{2} = [1 12 0 0; 2 0 12 0; -5 0 0 0]; 
    eqs{3} = [1 6 0 6; -1 0 0 0]; 

    system = systemstruct(eqs); 
end