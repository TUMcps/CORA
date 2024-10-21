function res = test_zonotope_enclosePoints
% test_zonotope_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = test_zonotope_enclosePoints
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% points
p = [1  3 -2  4 3 -1 1 0;...
     2 -1  1 -3 2  1 0 1];

% compute enclosing zonotope
Z = zonotope.enclosePoints(p);
% different method
Z_ = zonotope.enclosePoints(p,'stursberg');

% visualization
% figure; hold on;
% plot(Z);
% plot(p(1,:),p(2,:),'.k');

% check if all points are contained
assert(all(contains(Z,p)))
assert(all(contains(Z_,p)))

% ------------------------------ END OF CODE ------------------------------
