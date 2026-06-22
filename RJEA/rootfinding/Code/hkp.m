function [ma] = hkp(mep)
    % HKP computes the Hochstenbach-Kosir-Plestenjak (HKP) bound of a
    % multiparameter eigenvalue problem (MEP).
    %
    %   [ma] = HKP(mep) computes the HKP bound ma of mep.
    %
    %   Input/output variables:
    %  
    %       ma: [int] HKP bound.
    %       mep: [struct] multiparameter eigenvalue problem.
    
    % MacaulayLab (2022) - Christof Vermeersch.

    % TODO: consider the bound when k > l + n - 1.

    [~,l] = size(mep.P{1});
    ma = mep.dmax^mep.n*nchoosek(mep.n+l-1,mep.n);
end