function E = squareExpMat(E)
% squareExpMat - compute the exponent matrix for the squared polynomial zonotope
%
% Syntax:
%    E = nnHelper.squareExpMat(E)
%
% Inputs:
%    E - exponential matrix
%
% Outputs:
%    E - squared exponential matrix
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

E_ = [];

for i = 1:size(E, 2) - 1
    for j = i + 1:size(E, 2)
        E_ = [E_, E(:, i) + E(:, j)];
    end
end

E = [2 * E, E_];
end

%------------- END OF CODE --------------