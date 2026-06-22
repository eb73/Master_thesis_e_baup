classdef mepstruct < problemstruct
    methods
        function obj = mepstruct(mat,dmax,n,basis)
            if ~exist('basis','var')
                basis = 'monomial';
            else
                basis = basis;
            end
            s = numel(mat);
            supp = monomialsmatrix(dmax,n);
            [nrows,ncols] = size(mat{1});
            coef = zeros(s,nrows,ncols);
            for k = 1:s
                coef(k,:,:) = mat{k};
            end
            obj@problemstruct({coef},{supp},basis)
            obj.type = 'multiparameter eigenvalue problem';
        end 
    end
end