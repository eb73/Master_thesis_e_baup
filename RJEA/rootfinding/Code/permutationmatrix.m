function [P] = permutationmatrix(d,n,rows,shift,blocksize)
    % PERMUTATIONMATRIX constructs the permutation matrix used in the
    % column space based approach.
    %
    %   [P] = PERMUTATIONMATRIX(d,n,rows,shift) creates a permutation 
    %   matrix for a problem in which given rows are shifted by a certain 
    %   shift polynomial.
    %
    %   [___] = PERMUTATIONMATRIX(___,blocksize) considers the block 
    %   version instead.
    %
    %   Input arguments:
    %       d [int]: degree of the multivariate (block) Vandermonde matrix.
    %       n [int]: number of variables/eigenvalues.
    %       rows [int(ma)]: rows to shift.
    %       shift [double(k,n+1)]: shift polynomial.
    %       blocksize [int - optional]: size of the eigenvector.
    %
    %   Output argument:
    %       P [double(q,q)]: permutation matrix.
    %
    %   See also BASISMON, MONOMIALSMATRIX.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    if ~exist('blocksize','var')
        blocksize = 1;
    end

    p = nbmonomials(d,n)*blocksize;
    P = zeros(p,p);

    % Create B matrix:
    P(1:length(rows),rows) = eye(length(rows));
    rep = rows;

    Sg = shiftmatrix(d,n,rows,shift,blocksize);
    mask = ones(size(Sg,2),1);
    mask(rows) = zeros(length(rows),1);
    B2 = [];
    for k = 1:length(rows)
        remainder = Sg(k,find(mask));
        if any(remainder)
            B2 = [B2 rows(k)];
        end
    end

    % Create C matrix:
    P(length(rows)+1:length(rows)+length(B2),:) = Sg(B2,:);
    for k = 1:length(B2)
        rowsinsg = find(Sg(B2(k),:));
        for l = 1:length(rowsinsg)
            if ~ismember(rowsinsg(l),rep)
                rep = [rep rowsinsg(l)];
                break
            end
        end
    end
    notrep = setdiff(1:p,rep);

    % Create D matrix:
    P(length(rows)+length(B2)+1:end,notrep) = eye(length(notrep));     
end