function [mep] = muhic4()
    A = [2 1 0; 0 2 0; 0 1 2];
    B = [1 2 1; 1 2 0; 0 0 0];
    C = [0 0 0; 0 2 1; 1 2 1];
    D = [1 1 0; 0 1 0; 0 1 1];
    E = [1 2 1; 1 2 0; 0 0 0];
    F = [0 0 0; 0 2 1; 1 2 1];
    
    mat = cell(3,1);
    mat{1} = [kron(A,F) - kron(C,D); kron(B,D) - kron(A,E)];
    mat{2} = [kron(C,E) - kron(B,F); zeros(9,9)];
    mat{3} = [zeros(9,9); kron(C,E) - kron(B,F)];

    mep = mepstruct(mat,1,2);
end