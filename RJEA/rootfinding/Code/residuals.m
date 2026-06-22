function [accuracy, residualvalues] = residuals(problem,X,evalimp)  
    % RESIDUALS computes the residuals of a problem.
    %
    %   [accuracy, residualvalues] = RESIDUALS(problem,X) computes the
    %   residual values and accuracy (maximum residual values) of a problem
    %   in the solutions
    %
    %   [___] = RESIDUALS(___,evalimp) uses a user-defined evaluation
    %   function.
    %
    %   Input arguments:
    %       problem [struct]: system of multivariate polynomial equations
    %       or multiparameter eigenvalue problem.
    %       X [double(ma,n)]: affine solutions.
    %       evalimp [function - optional]: evaluation of the problem
    %       (default = @evalmon).
    %
    %   Output arguments:
    %       accuracy [double]: maximum residual value.
    %       residualvalues [double(ma,1)]: residual value of every affine
    %       solution.
    %
    %   See also EVALMON, EVALCHEB.

    % MacaulayLab (2023) - Christof Vermeersch.

    % Set the default basis:
    if ~exist('evalimp','var') 
        evalimp = @evalmon;
    end
    
    % Compute the residual errors:
    ma = size(X,1);
    residualvalues = zeros(ma,1);
    
    for sol = 1:ma
        singularvalues = svd(squeeze(problemvalue(problem,X(sol,:),evalimp)));
        residualvalues(sol) = singularvalues(end);
    end
    
    % Determine the maximum residual error:
    accuracy = max(residualvalues);
end