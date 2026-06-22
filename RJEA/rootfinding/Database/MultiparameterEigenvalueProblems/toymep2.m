function [mep] = toymep2()
    mat = cell(6,1);
    mat{1} = [1 2; 3 4; 3 4];
    mat{2} = [2 1; 0 1; 1 3];
    mat{3} = zeros(3,2);
    mat{4} = zeros(3,2);
    mat{5} = [3 4; 2 1; 0 1];
    mat{6} = [1 2; 4 2; 2 1];
    mep = mepstruct(mat,2,2);
end