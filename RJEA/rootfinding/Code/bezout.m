function [mb] = bezout(system)
    % BEZOUT computes the Bezout number of a system of multivariate 
    % polynomial equations.
    %
    %   [mb] = BEZOUT(system) computes the Bezout number of a system of
    %   multivariate polynomial equations.
    %
    %   Input argument:
    %       system [struct]: system of multivariate polynomial equations.
    %
    %   Output argument:
    %       mb [int]: Bezout number.
    %
    %   See also KUSHNIRENKO, BKK.
    
    % MacaulayLab (2023) - Christof Vermeersch.
    
    % TODO: add checks for the requirements of the Bezout number.
    % TODO: throw an error when MEP is given as input.

    mb = prod(vertcat(system.di{:}));
end