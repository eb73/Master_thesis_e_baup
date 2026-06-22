function [system] = toy1_cheb()
    eqs = cell(2,1); 
    eqs{1} = [-1.5 0 0; 1 1 0; -1.5 0 2]; 
    eqs{2} = [2 1 0; -6 0 1]; 

    system = systemstruct(eqs);
end