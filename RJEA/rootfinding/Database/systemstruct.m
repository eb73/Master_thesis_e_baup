classdef systemstruct < problemstruct
    methods
        function obj = systemstruct(eqs,basis)
            if ~exist('basis','var')
                basis = 'monomial';
            else
                basis = basis;
            end
            s = numel(eqs);
            coef = cell(s,1);
            supp = cell(s,1);
            for k = 1:s
                polynomial = eqs{k};
                coef{k} = polynomial(:,1);
                supp{k} = polynomial(:,2:end);
            end
            obj@problemstruct(coef,supp,basis)
            obj.type = 'system';
        end 
    end
end