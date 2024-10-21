function example_manual_example_matZonotope()
% example_manual_example_matZonotope - example from the manual demonstrating
%   matrix zonotopes
%
% Syntax:
%   example_manual_example_matZonotope()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       07-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Zcenter = [1, 2; 3, 4]; % center of matrix zonotope Z1
Zdelta(:, :, 1) = [1, 0; 1, 1]; % generators of matrix zonotope Z1
matZ1 = matZonotope(Zcenter, Zdelta); % instantiate matrix zonotope Z1

Zcenter = [-1, 2; 2, -1]; % center of matrix zonotope Z2
Zdelta(:, :, 1) = [0, 0.5; 0.5, 0]; % generators of matrix zonotope Z2
matZ2 = matZonotope(Zcenter, Zdelta); % instantiate matrix zonotope Z2

matZ3 = matZ1 + matZ2 % perform Minkowski addition and display result
matZ4 = matZ1 * matZ2 % compute multiplication of and display result

intZ = intervalMatrix(matZ1) % compute interval matrix and display result

% ------------------------------ END OF CODE ------------------------------
