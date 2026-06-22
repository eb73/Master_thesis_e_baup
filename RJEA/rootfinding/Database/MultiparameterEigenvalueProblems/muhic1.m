function [mep] = muhic1()
    A = [0 -1; -1 1];
    B = [-1 0; 0 0];
    C = [-1 0; 0 0];
    D = [0 -1; -1 2];
    E = [-2 0; 0 0];
    F = [1 0; 0 0];
    
    mat = cell(3,1);
    mat{1} = [kron(A,F) - kron(C,D); kron(B,D) - kron(A,E)];
    mat{2} = [kron(C,E) - kron(B,F); zeros(9,9)];
    mat{3} = [zeros(9,9); kron(C,E) - kron(B,F)];

    mep = mepstruct(mat,1,2);
end