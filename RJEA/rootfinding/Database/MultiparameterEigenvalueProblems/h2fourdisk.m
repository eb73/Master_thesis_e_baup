function [mep] = h2fourdisk()
    b = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219];
    a = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513];

    mat = cell(15,1);
    mat{1} = [b zeros(1,4); a zeros(1,3); 0 a 0 0; 0 0 a 0; zeros(1,3) a; zeros(2,8) eye(2)]';
    mat{2} = [zeros(1,4) b; zeros(4,10); zeros(2,4) -2*eye(2) zeros(2,4)]';
    mat{3} = [zeros(1,3) b 0; zeros(4,10); zeros(2,5) -2*eye(2) zeros(2,3)]';
    mat{4} = [zeros(1,2) b zeros(1,2); zeros(4,10); zeros(2,6) -2*eye(2) zeros(2,2)]';
    mat{5} = [0 b zeros(1,3); zeros(4,10); zeros(2,7) -2*eye(2) zeros(2,1)]';
    mat{6} = [zeros(5,10); -1*eye(2) zeros(2,8)]';
    mat{7} = [zeros(5,10); zeros(2,1) -2*eye(2) zeros(2,7)]';
    mat{8} = [zeros(5,10); zeros(2,2) -2*eye(2) zeros(2,6)]';
    mat{9} = [zeros(5,10); zeros(2,3) -2*eye(2) zeros(2,5)]';
    mat{10} = [zeros(5,10); zeros(2,2) -1*eye(2) zeros(2,6)]';
    mat{11} = [zeros(5,10); zeros(2,3) -2*eye(2) zeros(2,5)]';
    mat{12} = [zeros(5,10); zeros(2,4) -2*eye(2) zeros(2,4)]';
    mat{13} = [zeros(5,10); zeros(2,4) -1*eye(2) zeros(2,4)]';
    mat{14} = [zeros(5,10); zeros(2,5) -2*eye(2) zeros(2,3)]';
    mat{15} = [zeros(5,10); zeros(2,6) -1*eye(2) zeros(2,2)]';
    mep = mepstruct(mat,2,4);
end