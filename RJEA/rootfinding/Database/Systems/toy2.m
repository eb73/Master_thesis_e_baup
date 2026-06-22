function system = toy2()
    eqs = cell(2,1);
    eqs{1} = [1 2 0; -3 0 2; 1 0 0];
    eqs{2} = [3 2 0; -1.5 0 2; -1 4 0; 1 0 4];
    system = systemstruct(eqs);
end