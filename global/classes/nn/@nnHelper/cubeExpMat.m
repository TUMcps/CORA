function E = cubeExpMat(E)
% cubeExpMat - compute the exponent matrix for the cubed polynomial zonotope
%
% Syntax:
%    E = nnHelper.cubeExpMat(E)
%
% Inputs:
%    E - exponential matrix
%
% Outputs:
%    E - cubed exponential matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

%------------- BEGIN CODE --------------

E1 = [];
E2 = [];
E3 = [];

for i = 1:size(E, 2) - 1
    for j = i + 1:size(E, 2)
        E1 = [E1, 2 * E(:, i) + E(:, j)];
    end
end

for i = 2:size(E, 2)
    for j = 1:i - 1
        E2 = [E2, 2 * E(:, i) + E(:, j)];
    end
end

for i = 1:size(E, 2) - 1
    for j = i + 1:size(E, 2) - 1
        for k = j + 1:size(E, 2)
            E3 = [E3, E(:, i) + E(:, j) + E(:, k)];
        end
    end
end

E = [3 * E, E1, E2, E3];
end

%------------- END OF CODE --------------