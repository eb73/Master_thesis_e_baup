function [system] = kotsireas3()
    eqs = cell(3,1);
    eqs{1} = [2 3 3 0; -2 3 0 0; -4 2 3 1; 5 1 3 3; -1 0 3 5];
    eqs{2} = [-2 3 0 1; -3 1 3 2; 2 1 3 0; 1 0 3 4; -2 0 3 2];
    eqs{3} = [-4 1 0 2; 4 1 0 0; 1 0 2 0; 1 0 0 4; -2 0 0 2; 1 0 0 0];

    system = systemstruct(eqs);
end