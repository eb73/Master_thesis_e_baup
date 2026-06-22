function our_solutions = fixschur(A)
    d = size(A,1);
    [Q,~] = schur(A{1},'complex');
    D = cell(d,1);
    D{1} = diag(Q' * A{1}* Q);
    our_solutions = zeros(size(A{1},1),d-1);
    for k = 2:d
        our_solutions(:,k-1) =diag(Q' * A{k}* Q);
        D{k} = our_solutions(:,k-1);
    end
    for k = 1:d-1
        our_solutions(:,k) = D{k+1};
    end
end