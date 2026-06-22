function [system] = bifurcation() 
    eqs = cell(3,1); 
    eqs{1} = [5 9 0 0; -6 5 1 0; 1 1 4 0; 2 1 0 1]; 
    eqs{2} = [-2 6 1 0; 2 2 3 0; 2 0 1 1]; 
    eqs{3} = [1 2 0 0; 1 0 2 0; -0.265625 0 0 0]; 

    system = systemstruct(eqs); 
end