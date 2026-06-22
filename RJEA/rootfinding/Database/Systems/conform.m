function [system] = conform() 
    eqs = cell(3,1); 
    eqs{1} = [-3 0 2 2; -1 0 2 0; 8 0 1 1; -1 0 0 2; -9 0 0 0]; 
    eqs{2} = [-3 2 0 2; -1 2 0 0; 8 1 0 1; -1 0 0 2; -9 0 0 0]; 
    eqs{3} = [-3 2 2 0; -1 2 0 0; 8 1 1 0; -1 0 2 0; -9 0 0 0]; 

    system = systemstruct(eqs); 
end