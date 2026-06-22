function [system] = dreesen12() 
    % DREESEN12 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN12() returns the system of multivariate 
    %   polynomial equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 3 0; 1 0 3; -9 2 1; 20 1 1; -3 1 0; -20 0 0]; 
    eqs{2} = [1 2 0; 4 0 2; -1 1 1; -80 0 0]; 

    system = systemstruct(eqs); 
end