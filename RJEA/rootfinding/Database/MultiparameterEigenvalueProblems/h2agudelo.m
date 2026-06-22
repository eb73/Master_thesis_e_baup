function [mep] = h2agudelo()
    mat = cell(3,1);
    mat{1} = [1 -1 -1 0; 9 -12 0 -1; -10 -49 0 0; 0 -78 0 0];
    mat{2} = [0 0 0 0; 1 0 2 0; 9 0 0 2; -10 0 0 0];
    mat{3} = [zeros(2,4); zeros(2,2); -1*eye(2)];

    mep = mepstruct(mat,2,1);
end