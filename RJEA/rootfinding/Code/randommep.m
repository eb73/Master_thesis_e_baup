function mep = randommep(dmax,n,k,l)
    % RANDOMMEP generates a dense multiparameter eigenvalue parameter
    % problem with random coefficient matrices.
    %
    %   [mep] = RANDOMMEP(dmax,n,k,l) generates a dense multiparameter 
    %   eigenvalue problem with n eigenvalues and random k x l coefficient
    %   matrices (maximum total degree dmax).
    %
    %   Input/output variables:
    %
    %       mep: [struct] multiparameter eigenvalue problem.
    %       n: [int] number of eigenvalues.
    %       k: [int] number of rows of coefficient matrices.
    %       l: [int] number of columns of coefficient matrices.
    %       dmax: [int] maximum total degree.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    s = nbmonomials(dmax,n);
    P = cell(s,1);
    for i = 1:s
        P{i} = randn(k,l);
    end
    mep = mepstruct(P,dmax,n);
end