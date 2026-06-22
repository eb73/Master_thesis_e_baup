function [mep] = alsubaie2()    
    mat = cell(3,1);
    mat{1} = [2 3; 2 4; 4 4; 4 2];
    mat{2} = [3 3; 1 2; 2 3; 3 4];
    mat{3} = [3 1; 3 3; 3 4; 3 2];
    mat{4} = [4 2; 2 1; 4 1; 2 2];
    mep = mepstruct(mat,1,3);
end