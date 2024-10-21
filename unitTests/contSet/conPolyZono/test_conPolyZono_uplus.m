function res = test_conPolyZono_uplus
% test_conPolyZono_uplus - unit test function of uplus
%
% Syntax:
%    res = test_conPolyZono_uplus
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
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
c = [0;0];
G = [2 0 2; 0 2 2];
GI = [1;6];
E = [1 0 1; 0 1 1; 0 0 0];
A = [2 2 4 -4];
b = 0;
EC = [1 0 1 0; 0 1 1 0; 0 0 0 1];
cPZ = conPolyZono(c,G,E,A,b,EC,GI);

pcPZ = +cPZ;

% compare with cPZ
assert(isequal(pcPZ, cPZ));

% test empty case
assert(isemptyobject(+conPolyZono.empty(2)));

% add results
res = true;

% ------------------------------ END OF CODE ------------------------------
