function [system] = emiris() 
    eqs = cell(3,1); 
    eqs{1} = [-9 0 0 0; -1 0 2 0; -1 0 0 2; -3 0 2 2; 8 0 1 1]; 
    eqs{2} = [-9 0 0 0; -1 0 0 2; -1 2 0 0; -3 2 0 2; 8 1 0 1]; 
    eqs{3} = [-9 0 0 0; -1 2 0 0; -1 0 2 0; -3 2 2 0; 8 1 1 0]; 

    system = systemstruct(eqs); 
end