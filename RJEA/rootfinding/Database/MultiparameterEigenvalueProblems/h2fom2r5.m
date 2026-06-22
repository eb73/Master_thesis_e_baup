function [mep] = h2fom2r5()
    b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
    a = [1, 10, 46, 130, 239, 280, 194, 60];

    mat = cell(21,1);
    mat{1} = [b zeros(1,6); a zeros(1,4); zeros(1,1) a zeros(1,3); zeros(1,2) a zeros(1,2); zeros(1,3) a zeros(1,1); zeros(1,4) a;   -1*eye(2) zeros(2,10)]';
    mat{2} = [zeros(1,5) b; zeros(5,12); zeros(2,5) 2*eye(2) zeros(2,5)]';
    mat{3} = [zeros(1,4) b zeros(1,1); zeros(5,12); zeros(2,4) -2*eye(2) zeros(2,6)]';
    mat{4} = [zeros(1,3) b zeros(1,2); zeros(5,12); zeros(2,3) 2*eye(2) zeros(2,7)]';
    mat{5} = [zeros(1,2) b zeros(1,3); zeros(5,12); zeros(2,2) -2*eye(2) zeros(2,8)]';
    mat{6} = [zeros(1,1) b zeros(1,4); zeros(5,12); zeros(2,1) 2*eye(2) zeros(2,9)]';
    mat{7} = [zeros(6,12); zeros(2,10) -1*eye(2)];
    mat{8} = [zeros(6,12); zeros(2,9) 2*eye(2) zeros(2,1)];
    mat{9} = [zeros(6,12); zeros(2,8) -2*eye(2) zeros(2,2)];
    mat{10} = [zeros(6,12); zeros(2,7) 2*eye(2) zeros(2,3)];
    mat{11} = [zeros(6,12); zeros(2,6) -2*eye(2) zeros(2,4)];
    mat{12} = [zeros(6,12); zeros(2,8) -1*eye(2) zeros(2,2)];
    mat{13} = [zeros(6,12); zeros(2,7) 2*eye(2) zeros(2,3)];
    mat{14} = [zeros(6,12); zeros(2,6) -2*eye(2) zeros(2,4)];
    mat{15} = [zeros(6,12); zeros(2,5) 2*eye(2) zeros(2,5)];
    mat{16} = [zeros(6,12); zeros(2,6) -1*eye(2) zeros(2,4)];
    mat{17} = [zeros(6,12); zeros(2,5) 2*eye(2) zeros(2,5)];
    mat{18} = [zeros(6,12); zeros(2,4) -2*eye(2) zeros(2,6)];
    mat{19} = [zeros(6,12); zeros(2,4) -1*eye(2) zeros(2,6)];
    mat{20} = [zeros(6,12); zeros(2,3) 2*eye(2) zeros(2,7)];
    mat{21} = [zeros(6,12); zeros(2,2) -1*eye(2) zeros(2,8)];

    mep = mepstruct(mat,2,5);
end