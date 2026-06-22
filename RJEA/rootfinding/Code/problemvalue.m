function [values] = problemvalue(problem,x,evalimp)
    % PROBLEMVALUE evaluates the problem in a certain point.
    %
    %   [values] = PROBLEMVALUE(problem,x) computes the value of a problem
    %   in a given point.
    %
    %   [___] = PROBLEMVALUE(___,evalimp) uses a user-defined evaluation
    %   function.
    %
    %   Input arguments:
    %       problem [struct]: system of multivariate polynomial equations
    %       or multiparameter eigenvalue problem.
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

    [~,p,q] = size(problem.coef{1});

    % Evaluate all the polynomials:
    values = zeros(problem.s,p,q);
    for k = 1:problem.s
        values(k,:,:) = value(problem.coef{k},problem.supp{k},x,evalimp);
    end
end