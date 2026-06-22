function [system] = dreesen11() 
    % DREESEN11 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN11() returns the system of multivariate 
    %   polynomial equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 2 0; 1 1 1; 4 0 0]; 
    eqs{2} = [2 0 3; 2 1 2; 8 0 0]; 

    system = systemstruct(eqs); 
end