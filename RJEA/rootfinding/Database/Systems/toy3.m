function system = toy3()
    eqs = cell(2,1);
    eqs{1} = [-2 0 0; 1 1 1; 1 2 0];
    eqs{2} = [-2 0 0; 1 1 1; 1 0 2];
    system = systemstruct(eqs);
end