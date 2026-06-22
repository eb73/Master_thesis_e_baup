function [system] = dreesen5() 
    % DREESEN5 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN5() returns the system of multivariate polynomial 
    %   equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 0 1; -1 2 0]; 
    eqs{2} = [1 0 1; -2 1 0]; 

    system = systemstruct(eqs); 
end