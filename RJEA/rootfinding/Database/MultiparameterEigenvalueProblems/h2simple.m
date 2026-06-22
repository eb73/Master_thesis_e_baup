function [mep] = h2simple()
    A = [-12 -49 -78; 1 0 0; 0 1 0];
    B = [1; 0; 0];
    C = [1 9 -10];

    mat = cell(6,1);
    mat{1} = [zeros(1,9); zeros(1,5) 2*C 0; zeros(1,4) -2 zeros(1,4); zeros(3,1) A zeros(3,1) -eye(3) zeros(3,1); zeros(1,8) -1; zeros(3,5) A -B];
    mat{2} = [zeros(2,9); -2 zeros(1,8); zeros(3,1) -eye(3) zeros(3,5); zeros(1,4) -2 zeros(1,4); zeros(3,5) -eye(3) zeros(3,1)];
    mat{3} = [0 2*C zeros(1,5); zeros(1,4) -2 zeros(1,4); zeros(8,9)];
    mat{4} = zeros(10,9);
    mat{5} = zeros(10,9);
    mat{6} = [-1 zeros(1,8); zeros(9,9)];
    mep = mepstruct(mat,2,2); 
end