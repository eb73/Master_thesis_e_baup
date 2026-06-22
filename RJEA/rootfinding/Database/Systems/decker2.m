function [system] = decker2() 
    eqs = cell(2,1); 
    eqs{1} = [1 1 0; 1 0 3]; 
    eqs{2} = [1 2 1; -1 0 4]; 

    system = systemstruct(eqs); 
end