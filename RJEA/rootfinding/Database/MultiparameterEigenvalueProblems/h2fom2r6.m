function [mep] = h2fom2r6()
    b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
    a = -[1, 10, 46, 130, 239, 280, 194, 60];

    mat = cell(28,1);
    mat{1} = [b zeros(1,6); a zeros(1,5); 0 a zeros(1,4); zeros(1,2) a zeros(1,3); zeros(1,3) a zeros(1,2); zeros(1,4) a 0; zeros(1,5) a; -1 zeros(1,12)]';
    mat{2} = [zeros(1,6) b; zeros(6,13); zeros(1,6) -2 zeros(1,6)]';
    mat{3} = [zeros(1,5) b zeros(1,1); zeros(6,13); zeros(1,5) 2 zeros(1,7)]';
    mat{4} = [zeros(1,4) b zeros(1,2); zeros(6,13); zeros(1,4) -2 zeros(1,8)]';
    mat{5} = [zeros(1,3) b zeros(1,3); zeros(6,13); zeros(1,3) 2 zeros(1,9)]';
    mat{6} = [zeros(1,2) b zeros(1,4); zeros(6,13); zeros(1,2) -2 zeros(1,10)]';
    mat{7} = [zeros(1,1) b zeros(1,5); zeros(6,13); zeros(1,1) 2 zeros(1,11)]';
    mat{8} = [zeros(7,13); zeros(1,12) -1]';
    mat{9} = [zeros(7,13); zeros(1,11) 2 0]';
    mat{10} = [zeros(7,13); zeros(1,10) -2 zeros(1,2)]';
    mat{11} = [zeros(7,13); zeros(1,9) 2 zeros(1,3)]';
    mat{12} = [zeros(7,13); zeros(1,8) -2 zeros(1,4)]';
    mat{13} = [zeros(7,13); zeros(1,7) 2 zeros(1,5)]';
    mat{14} = [zeros(7,13); zeros(1,10) -1 zeros(1,2)]';
    mat{15} = [zeros(7,13); zeros(1,9) 2 zeros(1,3)]';
    mat{16} = [zeros(7,13); zeros(1,8) -2 zeros(1,4)]';
    mat{17} = [zeros(7,13); zeros(1,7) 2 zeros(1,5)]';
    mat{18} = [zeros(7,13); zeros(1,6) -2 zeros(1,6)]';
    mat{19} = [zeros(7,13); zeros(1,8) -1 zeros(1,4)]';
    mat{20} = [zeros(7,13); zeros(1,7) 2 zeros(1,5)]';
    mat{21} = [zeros(7,13); zeros(1,6) -2 zeros(1,6)]';
    mat{22} = [zeros(7,13); zeros(1,5) 2 zeros(1,7)]';
    mat{23} = [zeros(7,13); zeros(1,6) -1 zeros(1,6)]';
    mat{24} = [zeros(7,13); zeros(1,5) 2 zeros(1,7)]';
    mat{25} = [zeros(7,13); zeros(1,4) -2 zeros(1,8)]';
    mat{26} = [zeros(7,13); zeros(1,4) -1 zeros(1,8)]';
    mat{27} = [zeros(7,13); zeros(1,3) 2 zeros(1,9)]';
    mat{28} = [zeros(7,13); zeros(1,2) -1 zeros(1,10)]';

    mep = mepstruct(mat,2,6);
end