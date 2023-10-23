function completed = example_intervalMatrix()
% example_intervalMatrix - example for interval matrices
%
% Syntax:
%    completed = example_intervalMatrix()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalMatrix

% Authors:       Niklas Kochdumper
% Written:       25-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

Mcenter = [1 2; 3 4]; % center of interval matrix M1
Mdelta = [1 0; 1 1]; % delta of interval matrix M1
intM1 = intervalMatrix(Mcenter, Mdelta); % instantiate interval matrix M1

Mcenter = [-1 2; 2 -1]; % center of interval matrix M2
Mdelta = [0 0.5; 0.5 0]; % delta of interval matrix M2
intM2 = intervalMatrix(Mcenter, Mdelta); % instantiate inerval matrix M2

intM3 = intM1 + intM2 % perform Minkowski addition and display result
intM4 = intM1 * intM2 % compute multiplication and display result

matZ = matZonotope(intM1) % compute matrix zonotope and display result

completed = 1;

% ------------------------------ END OF CODE ------------------------------
