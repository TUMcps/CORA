function res = testLong_spectraShadow_conPolyZono
% testLong_spectraShadow_conPolyZono - unit test function of
%    conPolyZono
%
% Syntax:
%    res = testLong_spectraShadow_conPolyZono
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

% Authors:       Adrian Kulmburg
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% init polytope
A = [2 1; -1 3; -2 -2; 1 -3];
b = ones(4,1);
P = polytope(A,b);

% convert to spectraShadow
SpS = spectraShadow(P);

% convert to polyZonotope
cPZ = conPolyZono(SpS);

% not much more can be done at this stage, since these two set
% representations are completely different
res = true;


% ------------------------------ END OF CODE ------------------------------
