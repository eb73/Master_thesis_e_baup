function [X,output] = solvemep(problem,dend,options)
    % SOLVEMEP is an interface to macaulaylab(). It is recommended to use
    % macaulaylab() directly.
    %
    % See also MACAULAYLAB.

    % MacaulayLab (2023) - Christof Vermeersch.

    warning('off','backtrace')
    warning('This method is deprecated, use macaulaylab() instead.')

    [X,output] = macaulaylab(problem,dend,options);
end



