function [mep] = hkp2()
    mat = cell(6,1);
    mat{1} = [1 2; 3 4; 3 1];
    mat{2} = [1 3; 5 1; 1 4];
    mat{3} = [4 1; 1 3; 4 1];
    mat{4} = [2 3; 1 1; 1 2];
    mat{5} = [1 1; 2 2; 3 3];
    mat{6} = [3 1; 3 2; 1 2];
    mep = mepstruct(mat,2,2);
end