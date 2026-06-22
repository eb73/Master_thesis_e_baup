function [system] = cyclohexane1() 
    eqs = cell(3,1); 
    eqs{1} = [1313 0 2 2; 959 0 2 0; 1389 0 1 1; 774 0 0 2; -310 0 0 0]; 
    eqs{2} = [1269 2 0 2; 917 2 0 0; 1451 1 0 1; 755 0 0 2; -365 0 0 0]; 
    eqs{3} = [1352 2 2 0; 837 2 0 0; 1655 1 1 0; 838 0 2 0; -413 0 0 0]; 

    system = systemstruct(eqs); 
end