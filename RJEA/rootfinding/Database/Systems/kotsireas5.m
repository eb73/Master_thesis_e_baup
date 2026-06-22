function [system] = kotsireas5()
    eqs = cell(3,1);
    eqs{1} = [4 3 3 0; -6 3 0 0; -12 2 3 1; 15 1 3 3; -3 1 3 1; -3 0 3 5; 1 0 3 3];
    eqs{2} = [-6 3 0 1; -9 1 3 2; 5 1 3 0; 3 0 3 4; -5 0 3 2];
    eqs{3} = [-12 1 0 2; 12 1 0 0; 4 0 2 0; 3 0 0 4; -6 0 0 2; 3 0 0 0];
    
    system = systemstruct(eqs);
end