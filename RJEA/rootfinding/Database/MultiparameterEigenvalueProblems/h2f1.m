function [mep] = h2f1()
    b = [-2.9239, -39.5525,-97.5270, -147.1508]; 
    a = [1, 11.9584, 43.9119, 73.6759, 44.3821];

    mat = cell(6,1);
    mat{1} = [b 0 0; a 0; 0 a; -1*eye(2) zeros(2,4)]';
    mat{2} = [0 0 b; zeros(2,6); zeros(2,2) -2*eye(2) zeros(2,2)]';
    mat{3} = [0 b 0; zeros(2,6); zeros(2,1) 2*eye(2) zeros(2,3)]';
    mat{4} = [zeros(3,6); zeros(2,4) -1*eye(2)]';
    mat{5} = [zeros(3,6); zeros(2,3) 2*eye(2) zeros(2,1)]';
    mat{6} = [zeros(3,6); zeros(2,2) -1*eye(2) zeros(2,2)]';

    mep = mepstruct(mat,2,2);
end