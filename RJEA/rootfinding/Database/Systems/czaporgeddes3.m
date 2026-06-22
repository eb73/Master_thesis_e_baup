function [system] = czaporgeddes3(a,b,c,d,e,f,g,h,k)
    eqs = cell(2,1);
    eqs{1} = [a 2 0; b 1 1; c 1 0; d 0 2; 2 0 1; f 0 0];
    eqs{2} = [b 2 0; 4*d 1 1; 2*e 1 0; g 0 2; h 0 1; k 0 0];

    system = systemstruct(eqs);
end