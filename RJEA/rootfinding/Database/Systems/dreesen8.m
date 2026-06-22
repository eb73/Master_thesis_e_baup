function [system] = dreesen8() 
    % DREESEN8 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN8() returns the system of multivariate polynomial 
    %   equations. 

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(3,1); 
    eqs{1} = [1 2 0 0; 5 1 1 0; 4 0 1 1; -10 0 0 0]; 
    eqs{2} = [1 0 3 0; 3 2 1 0; -12 0 0 0]; 
    eqs{3} = [1 0 0 3; 4 1 1 1; -8 0 0 0]; 

    system = systemstruct(eqs); 
end