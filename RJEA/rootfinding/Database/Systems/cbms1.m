function [system] = cbms1() 
    eqs = cell(3,1); 
    eqs{1} = [1 3 0 0; -1 0 1 1]; 
    eqs{2} = [-1 1 0 1; 1 0 3 0]; 
    eqs{3} = [-1 1 1 0; 1 0 0 3]; 

    system = systemstruct(eqs); 
end