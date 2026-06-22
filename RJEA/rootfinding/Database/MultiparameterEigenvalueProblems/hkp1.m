function [mep] = hkp1()
    mat = cell(3,1);
    mat{1} = [1 2; 3 4; 3 1];
    mat{2} = [1 3; 5 1; 1 4];
    mat{3} = [4 1; 1 3; 4 1];
    mep = mepstruct(mat,1,2); 
end