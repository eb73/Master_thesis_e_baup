function [mep] = h2spanos3r4()
    b = [0, 0, 0, 0.00001, 0.0110, 1];
    a = -[1, 0.2220, 22.1242, 3.5445, 122.4433, 11.3231, 11.1100];

    mat = cell(15,1);
    mat{1} = [b 0 0 0 0; a 0 0 0; 0 a 0 0; 0 0 a 0; 0 0 0 a; -1*eye(2) zeros(2,8)]';
    mat{2} = [zeros(1,4) b; zeros(4,10); zeros(2,4) -2*eye(2) zeros(2,4)]';
    mat{3} = [zeros(1,3) b 0; zeros(4,10); zeros(2,3) 2*eye(2) zeros(2,5)]';
    mat{4} = [zeros(1,2) b 0 0; zeros(4,10); zeros(2,2) -2*eye(2) zeros(2,6)]';
    mat{5} = [0 b zeros(1,3); zeros(4,10); zeros(2,1) -2*eye(2) zeros(2,7)]';
    mat{6} = [zeros(5,10); zeros(2,8) -1*eye(2)]';
    mat{7} = [zeros(5,10); zeros(2,7) 2*eye(2) zeros(2,1)]';
    mat{8} = [zeros(5,10); zeros(2,6) -2*eye(2) zeros(2,2)]';
    mat{9} = [zeros(5,10); zeros(2,5) 2*eye(2) zeros(2,3)]';
    mat{10} = [zeros(5,10); zeros(2,6) -1*eye(2) zeros(2,2)]';
    mat{11} = [zeros(5,10); zeros(2,5) 2*eye(2) zeros(2,3)]';
    mat{12} = [zeros(5,10); zeros(2,4) -2*eye(2) zeros(2,4)]';
    mat{13} = [zeros(5,10); zeros(2,4) -1*eye(2) zeros(2,4)]';
    mat{14} = [zeros(5,10); zeros(2,3) 2*eye(2) zeros(2,5)]';
    mat{15} = [zeros(5,10); zeros(2,2) -1*eye(2) zeros(2,6)]';

    mep = mepstruct(mat,2,4);
end