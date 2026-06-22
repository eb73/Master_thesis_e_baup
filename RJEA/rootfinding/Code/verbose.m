function [] = verbose(output,varargin)
    % VERBOSE is an auxiliary function to display information about the
    % MacaulayLab solver.
    
    % MacaulayLab (2023) - Christof Vermeersch.


    linewidth = 75;
    linestr = pad('', linewidth,'_');
    hdrpad = pad('', floor((linewidth-11)/2) - 1);

    switch output

        case 'system'

            system = varargin{1};
            options = varargin{2};

            % Update required strings for verbose option:
            if options.blocked 
                blocked = 'block row-wise';
            else
                blocked = 'row-wise';
            end
            
            fprintf('\n')
            fprintf('<strong>%s MacaulayLab %s</strong> \n', hdrpad, hdrpad)
            fprintf('<strong>%s</strong> \n \n', linestr)
            fprintf(['The system has %d equations in %d variables (maximum' ...
                ' total degree is %d). \n \n'], system.s, system.n, system.dmax);
            fprintf('The selected polynomial basis is the %s basis \n \n', ...
                options.basis)
            fprintf(['The solver tackles the system via the null space ' ... 
            'of the Macaulay matrix: \n']);
            fprintf(['  *  Building a basis matrix of the null space ' ...
            '(%s - %s) \n'], options.recursive, blocked);
            fprintf(['\n \t | degree \t | nullity \t | increase \t |' ...
            'gap zone \t | \n'])
            fprintf(['\t |--------------------------------------------' ...
            '-------------------| \n'])

        case 'mep'

            mep = varargin{1};
            options = varargin{2};

            [s,k,l] = size(mep.coef{1});

            % Update required strings for verbose option:
            if options.blocked 
                blocked = 'block row-wise';
            else
                blocked = 'row-wise';
            end

        fprintf('\n')
        fprintf('<strong>%s MacaulayLab %s</strong> \n', hdrpad, hdrpad)
        fprintf('<strong>%s</strong> \n \n', linestr)
        fprintf(['The MEP is a %d-parameter eigenvalue problem' ...
            ' (maximum total degree is %d). \n'], mep.n, mep.dmax);
        fprintf(['The %d coefficient matrices are %d x %d matrices.' ...
            '\n \n'], s, k, l);
        fprintf('The selected polynomial basis is the %s basis \n \n', ...
            options.basis)
        fprintf(['The solver tackles the system via the null space ' ... 
            'of the block Macaulay matrix: \n']);
        fprintf(['  *  Building a basis matrix of the null space ' ...
            '(%s - %s) \n'], options.recursive, blocked);
        fprintf(['\n \t | degree \t | nullity \t | increase \t |' ...
            'gap zone \t | \n'])
        fprintf(['\t |--------------------------------------------' ...
            '-------------------| \n'])

        case 'initial'
            fprintf('\t | %d \t \t | %d \t \t | %s \t \t | NaN \t \t | \n', ...
            varargin{1}, varargin{2}, 'x')

        case 'iteration'
            fprintf(['\t | %d \t \t | %d \t \t | %d \t \t | %d \t \t' ...
                ' | \n'], varargin{1}, varargin{2}, varargin{2} - varargin{3}, varargin{4})

        case 'final' 
            time = varargin{6};
            options = varargin{7};
            output = varargin{8};

            fprintf(['\t | %d \t \t | %d \t \t | %d \t \t | %d \t \t' ...
                ' | \n'], varargin{1}, varargin{2}, varargin{2} - varargin{3}, varargin{4})
            fprintf('\n')
            fprintf('\t required degree = %d \n', varargin{1})
            fprintf('\t gap degree = %d \n', varargin{5})
            fprintf('\t gap zone size = %d \n', varargin{4})
            fprintf('\t total number of solutions = %d \n', output.mb)
            fprintf('\t number of affine solutions = %d \n', output.ma)
            fprintf('\t required time (enlarge) = %.4f seconds \n', time(1))
            fprintf('\t required time (check) = %.4f seconds \n', time(2))
            if options.null
                fprintf('  *  Performing a column compression: \n')
            else
                fprintf('  *  Executing the shifts in the column space: \n')
            end

        case 'compression'
            time = varargin{1};
            output = varargin{2};

            fprintf('\t from %d to %d solutions \n', output.mb, output.ma)
            fprintf('\t required time = %.4f seconds \n', time(3))
            fprintf('  *  Executing the shifts in the null space: \n')

        case 'shift'
            time = varargin{1};
            options = varargin{2};

            fprintf('\t required time = %.4f seconds \n', time(4))
            if options.clustering
                fprintf('  *  Clustering the solution candidates: \n')
            end

        case 'roots'
            time = varargin{1};
            options = varargin{2};
            output = varargin{3};

            if options.clustering
                fprintf('\t required time = %.4f seconds \n', time(5))
            end
            fprintf('\n')
            fprintf(['The solver results in %d affine solution candidates ' ...
                'in %.4f seconds. \n'], length(output.residuals), sum(time));
            if options.clustering && output.ma > 0
                fprintf('\t accuracy before clustering = %.4e \n', ...
                    output.accuracybeforeclustering);
                fprintf('\t accuracy after clustering = %.4e \n', ...
                    output.accuracy);
                fprintf('\t cluster tolerance = %.4e \n', options.clustertol);
            else
                fprintf('\t accuracy = %.4e \n', output.accuracy);
            end
            fprintf('\t rank tolerance = %.4e \n', options.tol);
            fprintf('<strong>%s</strong> \n \n', linestr)
            
        case 'eigentuples'

            time = varargin{1};
            options = varargin{2};
            output = varargin{3};
              
        fprintf('\n')
        fprintf(['The solver results in %d affine solution candidates ' ...
            'in %.4f seconds. \n'], length(output.residuals), sum(time));
        if options.clustering
            fprintf('\t accuracy before clustering = %.4e \n', ...
                output.accuracybeforeclustering);
            fprintf('\t accuracy after clustering = %.4e \n', ...
                output.accuracy);
            fprintf('\t cluster tolerance = %.4e \n', options.clustertol);
        else
            fprintf('\t accuracy = %.4e \n', output.accuracy);
        end
        fprintf('\t rank tolerance = %.4e \n', options.tol);
        fprintf('<strong>%s</strong> \n \n', linestr)

    end 
end