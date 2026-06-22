function [system] = ojika3() 
    eqs = cell(3,1); 
    eqs{1} = [1 1 0 0; 1 0 1 0; 1 0 0 1; -1 0 0 0]; 
    eqs{2} = [0.2 3 0 0; 0.5 0 2 0; 0.5 0 0 2; -1 0 0 1; 0.5 0 0 0]; 
    eqs{3} = [1 1 0 0; 1 0 1 0; 0.5 0 0 2; -0.5 0 0 0]; 

    system = systemstruct(eqs); 
end