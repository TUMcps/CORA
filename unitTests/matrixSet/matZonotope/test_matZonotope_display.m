function res = test_matZonotope_display
% test_matZonotope_display - unit test function for display (only check for
%    runtime errors)
% 
% Syntax:
%    res = test_matZonotope_display
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty matrix zonotope
matZ = matZonotope()

% scalar
C = 0;
G{1} = 1; G{2} = -2;
matZ = matZonotope(C,G)

% nx1 vector
C = [0; 1; 1];
G{1} = [1; -1; -2]; G{2} = [-2; 0; 1];
matZ = matZonotope(C,G)

% matrix
C = [0 2; 1 -1; 1 -2];
G{1} = [1 1; -1 0; -2 1]; G{2} = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G)

% ------------------------------ END OF CODE ------------------------------
