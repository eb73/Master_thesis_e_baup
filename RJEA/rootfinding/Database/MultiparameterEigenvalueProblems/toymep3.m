function [mep] = toymep3()
    mat = cell(4,1);
    mat{1} = [2 3; 2 5; 0 1; 1 1];
    mat{2} = [1 0; 0 1; 1 1; 2 1];
    mat{3} = [4 2; 2 3; 3 1; 3 1];
    mat{4} = [1 2; 1 4; 2 1; 4 2];
    mep = mepstruct(mat,1,3);
end