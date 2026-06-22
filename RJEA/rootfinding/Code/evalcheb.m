function [value] = evalcheb(coef,supp,x)
    % EVALCHEB evaluates the problem in the Chebyshev polynomial basis.
    %
    %   [value] = EVALCHEB(coef,supp,x) evaluates the problem given by its
    %   internal representation (coefficients and support) in the Chebyshev
    %   polynomial basis.
    %
    %   Input argument:
    %       coef [double(k,l,m)]: different coefficient/coefficient matrices.
    %       supp [cell(k,n)]: support of the problem.
    %       x [double(1,n)]: point to evaluate.
    %
    %   Output argument:
    %       value [double]: evaluation.
    %
    %   See also BASISCHEB.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    [~, p, q] = size(coef);
    supp = evalchebsupp(supp,x).*ones(1,p,q);
    value = squeeze(sum(coef.*supp,1));
end

function [suppvalue] = evalchebsupp(supp,x)
    suppvalue = zeros(size(supp,1),1);
    for k = 1:size(supp,1)
        suppvalue(k) = transformsupp(supp(k,:),x);
    end
end

function [monomialvalue] = transformsupp(monomial,x)

    dmax = max(monomial);
    T = cell(dmax+1,1);
    T{1} = [1 0];
    T{2} = [1 1];
    for k = 2:dmax
        T1 = T{k}; T1(:,1) = 2*T1(:,1); T1(:,2) = T1(:,2) + 1;
        T2 = T{k-1}; T2(:,1) = -T2(:,1);
        T{k+1} = flatpoly([T1; T2]);
    end
    
    term = 1;
    for l = 1:length(monomial)
        chebmon = T{monomial(l)+1};
        term = term*evalmon(chebmon(:,1),chebmon(:,2),x(l));
    end
    monomialvalue = term;
end