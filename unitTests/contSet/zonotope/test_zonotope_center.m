function res = test_zonotope_center
% test_zonotope_center - unit test function of center
%
% Syntax:
%    res = test_zonotope_center
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness of test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty zonotope
Z = zonotope.empty(2);
c = center(Z);
assert(isempty(c) && isnumeric(c) && all(size(c) == [2, 0]));

% 2D, create zonotope
c = [1; 5]; G = [2 3 4; 6 7 8];
Z = zonotope(c,G);
assert(all(withinTol(c,center(Z))));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
