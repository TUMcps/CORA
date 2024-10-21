function res = test_conPolyZono_copy
% test_conPolyZono_copy - unit test function of copy
%
% Syntax:
%    res = test_conPolyZono_copy
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D conPolyZono
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);
cPZ_copy = copy(cPZ);
assert(isequal(cPZ,cPZ_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
