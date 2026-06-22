function [system] = dreesen9() 
    % DREESEN9 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN9() returns the system of multivariate polynomial 
    %   equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(3,1); 
    eqs{1} = [1 1 1 0; -3 0 0 0]; 
    eqs{2} = [1 2 0 0; -1 0 0 2; 1 1 0 1; -5 0 0 0]; 
    eqs{3} = [1 0 0 3; -2 1 1 0; 7 0 0 0]; 

    system = systemstruct(eqs); 
end