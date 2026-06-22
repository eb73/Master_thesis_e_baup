function [mep] = realization(y,n)

    % a = (1, alpha_1, alpha_2, ..., alpha_n)
    a = sym('a',[1, n]); % Symbolic coefficients
    a = [1,a];
    N = length(y);

        % Compose T coefficient matrix (SISO case!!)
    T = sym('T', [N-n, N]);
    for i = 1:N-n
        T(i,1:i-1) = zeros(1,i-1);
        T(i,i:i-1+length(a)) = flip(a);
        T(i,i+length(a):end) = zeros(1,N-(i-1+length(a)));
    end
    D = T*T.';

        % Compose system of eqs from TK85, (13):
    f = sym('f',[(n+1)*(N-n)+1,1]); % f = (f, f_alpha_1, ..., f_alpha_n, -1)^T
    f(end) = -1;
    F(1:N-n) = D*f(1:N-n)+T*y*f(end);
    for i = 1:n
        F(i*(N-n)+1:(i+1)*(N-n)) = diff(D,a(i+1))*f(1:N-n)+D*f(i*(N-n)+1:(i+1)*(N-n))...
            + diff(T,a(i+1))*y*f(end);
    end
    for i = 1:n
        F((n+1)*(N-n)+i) = y.'*diff(T,a(i+1)).'*f(1:N-n)+y.'*T.'*f(i*(N-n)+1:(i+1)*(N-n));
    end

        % Compose MEP:
    [mat] = translate(F.');

    mep = mepstruct(mat,2,n);
end

function [mat] = translate(P)

    npol = length(P); %Number of polynomial equations
    Var = symvar(P);
    nvar=length(Var);
    Var_str = char(Var);
    
    for k=1:npol
       t=feval(symengine,'poly2list',vpa(P(k)),Var_str);
       n_mon = length(t); %number of monomials
        ma = zeros(n_mon,nvar+1);
       for j=1:n_mon
           monomial = t(j);
           ma(j,1) = monomial(1);
           ma(j,2:end) = monomial(2);
       end 
       Pol{k} = ma;
    end
    
    % Find max degree of each variable:
    d_max = zeros(nvar,1);
    for poly_index = 1:npol
        poly = Pol{poly_index};
        for mon_index = 1:size(poly,1)
            for var_index = 1:nvar
                d_max(var_index) = max(d_max(var_index), poly(mon_index,1+var_index));
            end
        end
    end
    
    % Select linear variables and put them in eigenvector:
    EigV_indices = find(d_max == 1); Var_mep_indices = find(d_max ~= 1);
    EigV = Var(EigV_indices); Var_mep = setdiff(Var, EigV);
    
    % Create monomials (both for MEP-matrices and eigenvector):
    monomials_MEP = monomialsmatrix(max(d_max),length(Var_mep));
    monomials_EigV = monomialsmatrix(1,length(EigV));
    
    % Compose matrices of MEP one by one: 
    mat = cell(1,size(monomials_MEP,1));
    for matrix_index = 1:length(monomials_MEP(:,1))
        M = zeros(npol,length(EigV)+1);
     
        for poly_index = 1:npol
            pol = Pol{poly_index}; 
            
            for mon_index = 1:size(pol,1)
                % Check whether this is the correct matrix for this monomial:
                tmp = pol(mon_index,Var_mep_indices+1);
                if monomials_MEP(matrix_index,:) == tmp
                    % Find correct column of matrix: 
                    mon_eig = pol(mon_index, EigV_indices+1);
                    [tf_2, index_2] = ismember(mon_eig,monomials_EigV,'rows');
                    if ~tf_2 % Should be always true and unique
                        warning('Something went wrong')
                        return
                    else
                       % Fill entry in matrix
                       M(poly_index,index_2) = double(pol(mon_index,1));
                    end
                end
            end
        end
        
        mat{matrix_index} = M;
    end
end


