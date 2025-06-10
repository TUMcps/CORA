function res = test_conPolyZono_plus
% test_conPolyZono_plus - unit test function of plus
%
% Syntax:
%    res = test_conPolyZono_plus
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
% Written:       24-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init sets ---

c = [0;0];
G = [2 2;2 -1];
E = [1 0; 0 1];
A = [1 1];
b = 0;
EC = [2 0; 0 1];
cPZ = conPolyZono(c,G,E,A,b,EC);

E = ellipsoid(0.1*[2 1;1 2]);

% compute Minkowski sums
% (only very simple tests here, sophisticated tests are in testLong)
S = cPZ + cPZ; 
S = cPZ + E; 

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
