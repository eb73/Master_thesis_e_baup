function [mep] = alsubaie1()
    mat = cell(3,1);
    mat{1} = [2 3; 2 4; 4 4];
    mat{2} = [3 3; 1 2; 2 3];
    mat{3} = [3 1; 3 3; 3 4];
    mep = mepstruct(mat,1,2);
end