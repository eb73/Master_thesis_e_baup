function [system] = dreesen13() 
    % DREESEN13 contains a system of multivariate polynomial equations. 
    % 
    %   [system] = DREESEN13() returns the system of multivariate 
    %   polynomial equations. 
    %
    %   Note: This system of multivariate polynomial equations is a noisy
    %   over-constrained problem.

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(4,1); 
    eqs{1} = [-0.671429 0 0; -0.100961 1 0; 0.667031 1 1; 0.031771 3 0; ...
        -0.0299538 2 1; 0.034994 0 3]; 
    eqs{2} = [-0.6714 0 0; -0.100145 1 0; 0.668768 1 1; 0.035137 3 0; ...
        -0.30288 2 1; 0.036594 0 3]; 
    eqs{3} = [-1.001697 0 0; 0.013063 2 0; -0.015508 1 1; 0.048478 0 2]; 
    eqs{4} = [-1.000747 0 0; 0.01558 2 0; -0.012807 1 1; 0.050581 0 2]; 

    system = systemstruct(eqs); 
end