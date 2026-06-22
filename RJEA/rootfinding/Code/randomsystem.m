function system = randomsystem(s,dmax,n)
    % RANDOMSYSTEM generates a dense system with random coefficients.
    %
    %   [system] = RANDOMSYSTEM(s,dmax,n) generates a dense system with n
    %   variables and s multivariate polynomial equations (maximum total
    %   degree dmax).
    %
    %   Input/output variables:
    %
    %       system: [struct] system of multivariate polynomial equations.
    %       n: [int] number of variables.
    %       s: [int] number of equations.
    %       dmax: [int] maximum total degree.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    K = monomialsmatrix(dmax,n);
    P = cell(s,1);
    for k = 1:s
        P{k} = [randn(size(K,1),1) K];
    end
    system = systemstruct(P);
end