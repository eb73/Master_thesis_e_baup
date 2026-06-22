function [system] = kotsireas8(m)
    eqs = cell(3,1);
    eqs{1} = [2*m 3 3 0; -2 3 0 0; -4 2 3 1; 5 1 3 3; 3*m-3 1 3 1; -1 0 3 5; -m+1 0 3 3];
    eqs{2} = [2 3 0 1; 3 1 3 2; -m-1 1 3 0; -1 0 3 4; m+1 0 3 2];
    eqs{3} = [-4 1 0 2; 4 1 0 0; 1 0 2 0; 1 0 0 4; -2 0 0 2; 1 0 0 0];

    system = systemstruct(eqs);
end