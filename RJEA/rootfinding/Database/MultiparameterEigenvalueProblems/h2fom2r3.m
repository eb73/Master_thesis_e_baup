function [mep] = h2fom2r3()
    b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
    a = -[1, 10, 46, 130, 239, 280, 194, 60];

    mat = cell(10,1);
    mat{1} = [b 0 0 0; a 0 0; 0 a 0; 0 0 a; -1*eye(4) zeros(4,6)]';
    mat{2} = [zeros(1,3) b; zeros(3,10); zeros(4,3) 2*eye(4) zeros(4,3)]';
    mat{3} = [zeros(1,2) b zeros(1,1); zeros(3,10); zeros(4,2) -2*eye(4) zeros(4,4)]';
    mat{4} = [zeros(1,1) b zeros(1,2); zeros(3,10); zeros(4,1) 2*eye(4) zeros(4,5)]';
    mat{5} = [zeros(4,10); zeros(4,6) -1*eye(4)]';
    mat{6} = [zeros(4,10); zeros(4,5) 2*eye(4) zeros(4,1)]';
    mat{7} = [zeros(4,10); zeros(4,4) -2*eye(4) zeros(4,2)]';
    mat{8} = [zeros(4,10); zeros(4,4) -1*eye(4) zeros(4,2)]';
    mat{9} = [zeros(4,10); zeros(4,3) 2*eye(4) zeros(4,3)]';
    mat{10} = [zeros(4,10); zeros(4,2) -1*eye(4) zeros(4,4)]';

    mep = mepstruct(mat,2,3);
end