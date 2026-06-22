function [mep] = h2f6()
    b = [41, 50, 140];
    a = conv([1, 1, 1],[1, 10, 100]);

    mat = cell(6,1);
    mat{1} = [0 b 0 0; a 0; 0 a; -1*eye(2) zeros(2,4)]';
    mat{2} = [zeros(1,3) b; zeros(2,6); zeros(2,2) -2*eye(2) zeros(2,2)]';
    mat{3} = [zeros(1,2) b 0; zeros(2,6); zeros(2,1) 2*eye(2) zeros(2,3)]';
    mat{4} = [zeros(3,6); zeros(2,4) -1*eye(2)]';
    mat{5} = [zeros(3,6); zeros(2,3) 2*eye(2) zeros(2,1)]';
    mat{6} = [zeros(3,6); zeros(2,2) -1*eye(2) zeros(2,2)]';

    mep = mepstruct(mat,2,2);
end