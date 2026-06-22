function [system] = dreesen10() 
    % DREESEN10 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN10() returns the system of multivariate 
    %   polynomial equations. 
    %
    %   Note: This system of multivariate polynomial equations has a
    %   positive-dimensional solution set at infinity.

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(4,1); 
    eqs{1} = [1 1 0 0 0; 1 0 1 0 0; -1 0 0 0 0]; 
    eqs{2} = [1 1 0 1 0; 1 0 1 0 1]; 
    eqs{3} = [1 1 0 2 0; 1 0 1 0 2; -0.6666666667 0 0 0 0]; 
    eqs{4} = [1 1 0 3 0; 1 0 1 0 3]; 

    system = systemstruct(eqs); 
end