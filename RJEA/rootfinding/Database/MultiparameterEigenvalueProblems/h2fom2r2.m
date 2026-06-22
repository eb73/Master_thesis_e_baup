function [mep] = h2fom2r2()
    b = [0, 0, 0, 0.00001, 0.0110, 1];
    a = [1, 0.2220, 22.1242, 3.5445, 122.4433, 11.3231, 11.1100];


    mat = cell(6,1);
    mat{1} = [b 0 0; a 0; 0 a; -1*eye(4) zeros(4,4)]';
    mat{2} = [zeros(1,2) b; zeros(2,8); zeros(4,2) -2*eye(4) zeros(4,2)]';
    mat{3} = [zeros(1,1) b 0; zeros(2,8); zeros(4,1) 2*eye(4) zeros(4,3)]';
    mat{4} = [zeros(3,8); zeros(4,4) -1*eye(4)]';
    mat{5} = [zeros(3,8); zeros(4,3) 2*eye(4) zeros(4,1)]';
    mat{6} = [zeros(3,8); zeros(4,2) -1*eye(4) zeros(4,2)]';

    mep = mepstruct(mat,2,2);
end