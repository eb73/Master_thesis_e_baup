function [system] = czapsorgeddes4(a,b,c,d,e,f,g,h,k)
    eqs = cell(3,1);
    eqs{1} = [1 2 0 0; a 0 1 1; d 1 0 0; g 0 0 0];
    eqs{2} = [1 0 2 0; b 1 0 1; e 0 1 0; h 0 0 0];
    eqs{3} = [1 0 0 2; c 1 1 0; f 0 0 1; k 0 0 0];

    system = systemstruct(eqs);
end