classdef problemstruct

    properties
        coef
        supp
        n
        s
        di
        dmax
        nnze
        basis
        ma
        mb
        macaulaylabtime
        citationkey
        type
    end

    methods
        function obj = problemstruct(coef,supp,basis)
            if ~exist('basis','var')
                obj.basis = 'monomial';
            else
                obj.basis = basis;
            end
            obj.coef = coef;
            obj.supp = supp;
            obj.n = size(supp{1},2);
            obj.s = numel(coef);
            obj.di = cell(obj.s,1);
            obj.nnze = cell(obj.s,1);
            for k = 1:obj.s
                obj.di{k} = max(sum(supp{k},2));
                obj.nnze{k} = nnz(coef{k});
            end   
            obj.dmax = max(horzcat(obj.di{:}));
        end 
    end
end