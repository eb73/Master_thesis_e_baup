function [ma] = bkk(system)
    % BKK computes the Bernstein–Khovanskii–Kushnirenko (BKK) bound of a 
    % sytem of multivariate polynomial equations
    %
    %   [ma] = BKK(system) computes the BKK bound on the number of affine
    %   solutions of a system of multivariate polynomial equations.
    %
    %   Input argument:
    %       system [struct]: system of multivariate polynomial equations.
    %
    %   Output argument:
    %       ma [int]: BKK bound.
    %
    %   See also BEZOUT, KUSHNIRENKO.
    
    % MacaulayLab (2023) - Christof Vermeersch.

    % Note: this is a rather naive implementation.

    % TODO: add checks for the requirements of the BKK bound.
    % TODO: throw an error when MEP is given as input.
    % TODO: improve this naive implementation.

    % Compute the BKK bound, using Theorem 4.12 of Cox et al. (2007):
    ma = 0;
    for k = 1:system.n
        combinations = nchoosek(1:system.n,k);
        volume = 0;
        for l = 1:size(combinations,1)
            combo = combinations(l,:);
            B = cell(length(combo));
            for i = 1:length(combo)
                B{i} = system.supp{combo(i)};
            end
            [~, subvolume] = convhulln(minkowski(B));
            volume = volume + subvolume;
        end
        ma = ma + (-1)^(system.n-k)*volume;   
    end
end

function [sum] = minkowski(A)
    % TODO: determine requirements and implement this function correctly.
    sum = A{1};
    for k = 2:length(A)
        temp = [];
        for l = 1:size(A{k},1)
            temp = [temp; sum + A{k}(l,:)];
        end
        sum = unique(temp,'rows');
    end
end