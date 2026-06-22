function system = toy1()
    eqs = cell(2,1);
    eqs{1} = [7 0 0; -6 1 0; 1 2 0; 1 0 2];
    eqs{2} = [-3 0 0; 1 1 0; -1 0 1];
    system = systemstruct(eqs);
end