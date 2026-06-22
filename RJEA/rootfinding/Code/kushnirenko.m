function [ma] = kushnirenko(system)
    % KUSHNIRENKO computes the Kushnirenko bound of a system of 
    % multivariate polynomial equations.
    %
    %   [ma] = KUSHNIRENKO(system) computes the Kushnirenko bound on the
    %   number of affine solutions of a system of multivariate polynomial
    %   equations.
    %
    %   Input argument:
    %       system [struct]: system of multivariate polynomial equations.
    %
    %   Output argument:
    %       ma [int]: Kushnirenko bound.
    %       
    %   See also BEZOUT, BKK.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    % TODO: add checks for the requirements of the Kushnirenko bound.

    % Determine the unique support of the system:
    uniquesupport = unique(cat(1,system.supp{:}),'rows');

    % Compute the Kushnirenko bound:
    [~, volume] = convhulln(uniquesupport);
    ma = factorial(system.n)*volume;
end