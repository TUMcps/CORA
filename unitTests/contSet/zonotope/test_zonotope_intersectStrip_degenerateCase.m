function res = test_zonotope_intersectStrip_degenerateCase
% test_zonotope_intersectStrip_degenerateCase - unit test function of 
% intersectStrip to check whether a degenerate result is obtained.
% In a previous version, the method 'bravo' in [1] produced a result of 
% length 1e16.
%
% Syntax:
%    res = test_zonotope_intersectStrip_degenerateCase
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       23-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%% Simple 2D example which can be easily visualized
% zonotope
Z = zonotope([ ...
0.0070971205961507, 0.0100000000000000, -0.0200000000000000, -0.0000000000000000, 0.0340000000000000, 0.0420000000000000, 0.0160000000000000, -0.1200000000000000; ...
-0.0060883901920886, -0.0300000000000000, -0.0400000000000000, 0.0000000000000000, -0.1020000000000000, -0.1260000000000000, -0.0480000000000000, 0.0200000000000000]);


% strip
C = [-2, 1];
y = -0.1265718260823974;
phi = 0.2;

% obtain over-approximative zonotope after intersection
Z_over = intersectStrip(Z,C,phi,y,'bravo');

% check whether result is not too large
box = interval(Z_over);

% set enclosing box
box_encl = interval([-1;-1],[1;1]);

% check enclosure
res = contains(box_encl,box);

% ------------------------------ END OF CODE ------------------------------
