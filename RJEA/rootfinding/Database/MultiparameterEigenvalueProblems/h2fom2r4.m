function [mep] = h2fom2r4()
    b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
    a = -[1, 10, 46, 130, 239, 280, 194, 60];

    mat = cell(15,1);
    mat{1} = [b zeros(1,4); a zeros(1,3); 0 a zeros(1,2); zeros(1,2) a 0; zeros(1,3) a;  -1*eye(3) zeros(3,8)]';
    mat{2} = [zeros(1,4) b; zeros(4,11); zeros(3,4) -2*eye(3) zeros(3,4)]';
    mat{3} = [zeros(1,3) b zeros(1,1); zeros(4,11); zeros(3,3) -2*eye(3) zeros(3,5)]';
    mat{4} = [zeros(1,2) b zeros(1,2); zeros(4,11); zeros(3,2) -2*eye(3) zeros(3,6)]';
    mat{5} = [zeros(1,1) b zeros(1,3); zeros(4,11); zeros(3,1) -2*eye(3) zeros(3,7)]';
    mat{6} = [zeros(5,11); zeros(3,8) -1*eye(3)]';
    mat{7} = [zeros(5,11); zeros(3,7) 2*eye(3) zeros(3,1)]';
    mat{8} = [zeros(5,11); zeros(3,6) -2*eye(3) zeros(3,2)]';
    mat{9} = [zeros(5,11); zeros(3,5) 2*eye(3) zeros(3,3)]';
    mat{10} = [zeros(5,11); zeros(3,6) -1*eye(3) zeros(3,2)]';
    mat{11} = [zeros(5,11); zeros(3,5) 2*eye(3) zeros(3,3)]';
    mat{12} = [zeros(5,11); zeros(3,4) -2*eye(3) zeros(3,4)]';
    mat{13} = [zeros(5,11); zeros(3,4) -1*eye(3) zeros(3,4)]';
    mat{14} = [zeros(5,11); zeros(3,3) 2*eye(3) zeros(3,5)]';
    mat{15} = [zeros(5,11); zeros(3,2) -1*eye(3) zeros(3,6)]';

    mep = mepstruct(mat,2,4);
end