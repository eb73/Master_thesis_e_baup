function [s,n,di,dmax,nnze,k,l] = probleminfo(problem)
    % PROBLEMINFO retrieves the information about the problem.
    %
    %   [s,n,di,dmax,nnze,k,l] = PROBLEMINFO(problem) returns the basic
    %   properties of the problem.
    %
    %   Input argument:
    %       problem [struct]: system of multivariate polynomial equations
    %       or multiparameter eigenvalue problem.
    %
    %   Output arguments:
    %       s [int]: number of (matrix) equations.
    %       n [int]: number of variables.
    %       di [cell(s,1)]: total degree of every (matrix)
    %       equation.
    %       dmax [int]: maximum total degree of the problem.
    %       nnze [cell(s,1)]: number of non-zero coefficients of every
    %       (matrix) equation.
    %       k [int]: number of rows of the coefficient matrices.
    %       l [int]: number of columns of the coefficient matrices.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    s = problem.s;
    n = problem.n;
    di = problem.di;
    dmax = problem.dmax;
    nnze = problem.nnze;
    k = size(problem.coef{1},2);
    l = size(problem.coef{1},3);
end