function [system] = dreesen14() 
    % DREESEN14 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN14() returns the system of multivariate 
    %   polynomial equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 1 1; -8 0 0]; 
    eqs{2} = [1 2 0; -4 0 0]; 

    system = systemstruct(eqs); 
end