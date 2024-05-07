function res = test_matZonotope_expm
% test_matZonotope_expm - unit test function for expm
% 
% Syntax:
%    res = test_matZonotope_expm
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

% Authors:       Tobias Ladner
% Written:       26-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only test for runtime error
res = true;

% empty matrix zonotope
C = [0 1;0 -2.5];
D = [0 0;0 0.5];
intMat = intervalMatrix(C,D);
matZ = matZonotope(intMat);
eZ = expm(matZ);

eZ = expm(matZ,10);

% ------------------------------ END OF CODE ------------------------------
