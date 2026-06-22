function [mep] = toymep1()
    mat = cell(3,1);
    mat{1} = [2 6; 4 5; 0 1];
    mat{2} = [1 0; 0 1; 1 1];
    mat{3} = [4 2; 0 8; 1 1];
    mep = mepstruct(mat,1,2);
end