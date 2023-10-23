function res = test_conPolyZono_minkDiff
% test_conPolyZono_minkDiff - unit test function for the Minkowski
%    difference of constrained polynomial zonotopes
%
% Syntax:
%    res = test_conPolyZono_minkDiff()
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
% Written:       31-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% instantiate constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
GI = [0.5 -2; 1 1];
cPZ = conPolyZono(c,G,E,A,b,EC,GI);

% vector
shift = [-0.5;0];

% compute Minkowski difference (only center shifted)
cPZ_ = minkDiff(cPZ,shift);

if ~compareMatrices(cPZ_.c,c-shift) ...
        || ~compareMatrices([cPZ.G;cPZ.E],[cPZ_.G;cPZ_.E]) ...
        || ~compareMatrices([cPZ.A';cPZ.b']',[cPZ_.A';cPZ_.b']') ...
        || ~compareMatrices(cPZ.GI,cPZ_.GI)
    res = false;
end

% instantiate interval
I = interval([-0.02;-0.05],[0.04;0.03]);

% check if function runs through
cPZ_ = minkDiff(cPZ,I);

% ------------------------------ END OF CODE ------------------------------
