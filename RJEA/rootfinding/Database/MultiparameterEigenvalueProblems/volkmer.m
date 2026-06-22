function [mep] = volkmer()
    A = [-4 0 0; 0 0 0; 0 0 0];
    B = [1 0 0; 0 6 0; 0 0 1];
    C = [0 1 0; 1 0 1; 0 1 0];
    D = [-20 0; 0 0];
    E = [0 sqrt(3); sqrt(3) 0];
    F = [7 0; 0 1];
    
    mat = cell(3,1);
    mat{1} = [kron(A,F) - kron(C,D); kron(B,D) - kron(A,E)];
    mat{2} = [kron(C,E) - kron(B,F); zeros(9,9)];
    mat{3} = [zeros(9,9); kron(C,E) - kron(B,F)];

    mep = mepstruct(mat,1,2);
end