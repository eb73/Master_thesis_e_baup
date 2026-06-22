function [value] = value(coef,supp,x,evalimp)
    % VALUE evaluates the (matrix) equation in a certain point.
    %
    %   [value] = VALUE(coef,supp,x) computes the value of the (matrix)
    %   equation in a given point.
    %
    %   [___] = VALUE(___,evalimp) uses a user-defined evaluation function.
    %
    %   Input arguments:
    %       coef [double(k,l,m)]: different coefficient/coefficient matrices.
    %       supp [cell(k,n)]: support of the problem.
    %       x [double(1,n)]: point to evaluate.
    %       evalimp [function - optional]: evaluation of the problem
    %       (default = @evalmon).
    %
    %   Output arguments:
    %       values [double(s,k,l)]: evaluation of the problem.
    %
    %   See also EVALMON, EVALCHEB.

    % MacaulayLab (2023) - Christof Vermeersch.
    
    % Set the default basis:
    if ~exist('evalimp','var') 
        evalimp = @evalmon;
    end
    
    % Translate the affine point to a row vector:
    [nr,nc] = size(x);
    if nr > nc
        x = x.';
    end

    value = evalimp(coef,supp,x);
end