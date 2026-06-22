function [mep] = muhic2()
    A = [1 2 3 4 5 6; zeros(5,1) -eye(5)];
    B = [0 0 0 7 8 9; 1 zeros(1,5); zeros(1,6); 0 1 zeros(1,4); 0 0 1 zeros(1,3); zeros(1,6)];
    C = [0 0 0 0 0 10; zeros(1,6); 1 zeros(1,5); zeros(2,6); 0 0 1 zeros(1,3)];
    D = [10 9 8 7 6 5; zeros(5,1) -eye(5)];
    E = [0 0 0 4 3 2; 1 zeros(1,5); zeros(1,6); 0 1 zeros(1,4); 0 0 1 zeros(1,3); zeros(1,6)];
    F = [zeros(1,5) 1; zeros(1,6); 1 zeros(1,5); zeros(2,6); 0 0 1 zeros(1,3)];
    
    mat = cell(3,1);
    mat{1} = [kron(A,F) - kron(C,D); kron(B,D) - kron(A,E)];
    mat{2} = [kron(C,E) - kron(B,F); zeros(9,9)];
    mat{3} = [zeros(9,9); kron(C,E) - kron(B,F)];

    mep = mepstruct(mat,1,2);
end