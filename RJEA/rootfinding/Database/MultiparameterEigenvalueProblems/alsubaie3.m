function [mep] = alsubaie3()
    mat = cell(4,1);
    mat{1} = [2 2 3; 7 2 3; 2 2 9; 2 4 8; 1 4 6];
    mat{2} = [2 6 2; 3 4 3; 0 0 4; 0 0 1; 0 0 7];
    mat{3} = [5 4 7; 10 2 3; 0 0 6; 0 0 9; 0 0 6];
    mat{4} = [4 6 1; 3 8 3; 0 0 4; 0 0 7; 0 0 10];
    mep = mepstruct(mat,1,3);
end