function example_manual_example_intervalMatrix()
% example_manual_example_intervalMatrix - example from the manual demonstrating
%   interval matrix
%
% Syntax:
%   example_manual_example_intervalMatrix()
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

Mcenter = [1, 2; 3, 4]; % center of interval matrix M1
Mdelta = [1, 0; 1, 1]; % delta of interval matrix M1
intM1 = intervalMatrix(Mcenter, Mdelta); % instantiate interval matrix M1

Mcenter = [-1, 2; 2, -1]; % center of interval matrix M2
Mdelta = [0, 0.5; 0.5, 0]; % delta of interval matrix M2
intM2 = intervalMatrix(Mcenter, Mdelta); % instantiate interval matrix M2

intM3 = intM1 + intM2 % perform Minkowski addition and display result
intM4 = intM1 * intM2 % compute multiplication of and display result

matZ = matZonotope(intM1) % compute matrix zonotope and display result

% ------------------------------ END OF CODE ------------------------------
