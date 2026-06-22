function [X,output,commuting_family] = macaulaylab(problem,dend,options)
    % MACAULAYLAB solves a system of multivariate polynomial equations or
    % (rectangular) multiparameter eigenvalue problem via the (block)
    % Macaulay matrix.
    %
    %   [X] = MACAULAYLAB(problem) solves the problem with default options. 
    %   The final degree of the function is then set equal to 30.
    %   
    %   [___] = MACAULAYLAB(___,dend) uses a user-specified final degree.
    %   It is recommended to set this argument!
    %
    %   [___,output] = MACAULAYLAB(___,options) uses user-specified options
    %   and stores the output arguments.
    %
    %   Input arguments:
    %       problem [struct]: system of multivariate polynomial equations
    %       or (rectangular) multiparameter eigenvalue problem.
    %       dend [int - optional/recommended]: final degree (default = 30).
    %       options [struct - optional]: user-specified options.
    %           - options.null [boolean]: selects the null space based
    %             approach, instead of column space based approach (default
    %             = true).
    %           - options.blocked [boolean]: selects the blocked version of 
    %             the algorithm (default = true).
    %           - options.recursive [string]: selects an approach to build 
    %             the null space ('full', 'iterative', 'recursive', 
    %             default = 'sparse').
    %           - options.recrcheck [boolean]: selects the recursive rank
    %             checks (default = false).
    %           - options.tol [double]: tolerance in the rank checks 
    %             (default = 1e-10).
    %           - options.clustertol [double]: tolerance in the clustering 
    %             step (default = 1e-3).
    %           - options.shiftpoly [double(m,n+1)]: main shift polynomial 
    %             (default = [randn(n,1) eye(n)]).
    %           - options.verbose [boolean]: displays additional 
    %             information (default = false).
    %           - options.clustering [boolean]: choice to apply clustering 
    %             to gather the multiplicities (default = true).
    %           - options.posdim [boolean]: flag to solve problems with a 
    %             positive-dimensional solution set at infinity (default
    %             = false).
    %           - options.basis [string]: defines the polynomial basis
    %             (default = 'monomial', 'chebyshev', 'other').
    %           - options.mult [function]: contains the multiplication
    %             property of the 'other' polynomial basis.
    %           - options.eval [function]: contains an evaluation approach
    %             of a polynomial in the 'other' polynomial basis.
    %   
    %   Output arguments:
    %       X [double(ma,n)]: affine candidate solutions of the problem.
    %       output [struct]: additional output information.
    %           - output.accuracy [double]: contains the accuracy.
    %           - output.accuracybeforeclustering [double]: contains the 
    %             accuracy before the clustering step.
    %           - output.residuals [double(ma,1)]: contains the residual
    %             errors.
    %           - output.residualsbeforeclustering [double(ma,1)]: 
    %             contains the residual errors before the clustering step.
    %           - output.shiftvalues [double(ma,1)]: contains the 
    %             evaluations of the shift polynomial in all the affine 
    %             candidate solutions.
    %           - output.nullities [int(dend,1)]: contains the nullity of
    %             the Macaulay matrix for the computed degrees.
    %           - output.time [double(6,1)]: contains the required 
    %             computation time of the different part of the algorithm.
    %             In order to know the total required computation time,
    %             just use sum(output.time).
    %           - output.solutions [cell(n+1,1)]: contains the candidate
    %             solutions
    %       
    %   See also SHIFTNULLSPACE and SHIFTCOLUMNSPACE.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    % Fix the important variables:
    [s,n,~,dmax,~] = probleminfo(problem);
    blocksize = size(problem.coef{1},3);

    % Process the optional parameters:
    if ~exist('dend','var') 
        dend = 30; 
    end
    if ~exist('options','var') 
        options = struct();
    end
    
    parser = inputParser;
    parser.KeepUnmatched = true;
    parser.addParameter('blocked', true);
    parser.addParameter('recursive', 'sparse');
    parser.addParameter('shiftpoly', [randn(n,1) eye(n)]);
    parser.addParameter('verbose', false);
    parser.addParameter('tol', 1e-10);
    parser.addParameter('recrcheck', false);
    parser.addParameter('posdim', false);
    parser.addParameter('null', true);
    parser.addParameter('basis', 'monomial');
    parser.addParameter('clustering', false); 
    parser.addParameter('clustertol', 1e-3);
    parser.parse(options);
    options = parser.Results;

    % Display the properties of the system (verbose only):
    if options.verbose
        if blocksize == 1
            verbose('system',problem,options)
        else
            verbose('mep',problem,options)
        end
    end

    % Process the choice of polynomial basis:
    switch options.basis
        case 'monomial'
            options.mult = @basismon;
            options.eval = @evalmon;
        case 'chebyshev'
            options.mult = @basischeb;
            options.eval = @evalcheb;
    end
    
    % Determine the required gap degree:
    G = options.shiftpoly;
    dg = max(sum(G(:,2:end),2));
    dinit = max(dmax,dg);
    
    % Set the tolerance variable and stabilization flag:
    tol = options.tol;

    time = zeros(6,1);
    nullities = zeros(dend,1);
    gapsize = NaN;

    % Build the initial subspace:
    tic
    M = full(macaulay(problem,dinit,options.mult));
    if options.null
        Z = null(M);
        nullities(dinit) = size(Z,2);
    else
        Z = NaN;
        r = rank(M);
        nullities(dinit) = size(M,2) - r;
    end
    time(1) = toc;

    if options.verbose
        verbose('initial',dinit,nullities(dinit))
    end

    for d = dinit+1:dend
        % Enlarge the subspace:
        tic
        [M,Z] = enlarge(M,Z,problem,d,options);
        time(1) = time(1) + toc;

        % Check the subspace:
        tic
        if options.null
            nullities(d) = size(Z,2);
        else
            r = rank(full(M));
            nullities(d) = size(M,2) - r;
        end
        time(2) = time(2) + toc;

        if options.posdim || nullities(d) == nullities(d-dg)
            % Look for a gap zone:
            tic 
            [dgap,ma,gapsize,s,sinf] = check(M,Z,d,n,tol,blocksize,options);
            time(2) = time(2) + toc;

            % Break the iterative procedure when there is a gap zone:
            if gapsize >= dg
                if options.posdim
                    mb = Inf;
                else
                    mb = nullities(d);
                end
                output.mb = mb;
                output.ma = ma;
                output.size = size(M);
                break
            end
        end

        if options.verbose
            verbose('iteration',d,nullities(d),nullities(d-1),gapsize)
        end
        if d == dend
            error('No gap has been found for this final degree!')
        end
    end
    if options.verbose
        verbose('final',d,nullities(d),nullities(d-1),gapsize,dgap,time,options,output)
    end

    % Perform a column compression (only for null space based approach):
    if options.null && mb > ma
        tic
        Z = columncompression(Z,dgap + dg - 1,n,ma,blocksize);
        time(3) = toc;
        if options.verbose
            verbose('compression',time,output)
        end
    end

    % Solve the problem:
    tic
    if options.null
        [D, commuting_family] = shiftnullspace(Z,G,dgap,s,blocksize,options.mult);
    else 
        D = shiftcolumnspace(M,G,dgap,s,sinf,blocksize,options.mult);
    end
    time(4) = toc;
    if options.verbose
        verbose('shift',time,options)
    end

    % Check the accuracy of the obtained affine solutions:
    X = zeros(ma,n);
    for k = 1:n
        X(:,k) = D{k+1};
    end

    tic
    [accuracy, residualvalues] = residuals(problem,X,options.eval);
    time(6) = toc;

    % Perform the optional clustering.
    if options.clustering && ma > 0
        tic
        D = cluster(D,options.clustertol);
        time(5) = toc;
        output.accuracybeforeclustering = accuracy;
        output.residualvaluessbeforeclustering = residualvalues;
        for k = 1:n
            X(:,k) = D{k+1};
        end
        tic
        [accuracy, residualvalues] = residuals(problem,X,@options.eval);
        time(6) = time(6) + toc;
    end

    % Store the outputs:
    output.accuracy = accuracy;
    output.residuals = residualvalues;
    output.nullities = nullities;
    output.time = time;
    output.solutions = D;

    % Display the properties of the solution candidates (verbose only):
    if options.verbose
        if blocksize == 1
            verbose('roots',time,options,output)
        else
            verbose('eigentuples',time,options,output)
        end
    end
end

%% Define auxiliary functions:
function [M,Z] = enlarge(M,Z,problem,d,options)
    if options.null
        switch options.recursive
            case 'iterative'
                M = full(macaulay(problem,d,options.mult));
                Z = null(M);
            case 'recursive'
                rows = size(M,1);
                M = macaulayupdate(M,problem,d,options.mult);
                Z = nullrecrmacaulay(Z,full(M(rows+1:end,:)));
            otherwise
                M = NaN;
                Z = nullsparsemacaulay(Z,problem,d,options.tol,options.mult);
        end
    else
        switch options.recursive
            case 'iterative'
                M = full(macaulay(problem,d,options.mult));
            otherwise
                M = macaulayupdate(M,problem,d,options.mult);
        end
        Z = NaN;
    end
end

function [dgap,ma,gapsize,s,sinf] = check(M,Z,d,n,tol,vectorsize,options)
    if options.blocked && options.null % Find gap via blocked or non-blocked approach.
            % Set the standard monomials to NaN:
            s = NaN;
            sinf = NaN;

            [dgap,ma,gapsize] = gapzone(Z,d,n,tol,options.recrcheck,vectorsize);
    else
        if options.null
            % Find the standard monomials:
            s = stdmonomials(Z,tol);
        else
            % Find the standard monomials:
            sinv = stdmonomials(full(flip(M')));
            s = setdiff(1:size(M,2),size(M,2) - sinv + 1);
        end
        % Determine the number of affine solutions:
        [dgap,ma,gapsize] = gapzonestdmonomials(s,d,n,vectorsize);
        sinf = s(ma+1:end);
        s = s(1:ma);
    end
end



