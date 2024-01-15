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

res = true(0);

% empty zonotope
Z = zonotope.empty(2);
c = center(Z);
res(end+1,1) = isempty(c) && isnumeric(c) && all(size(c) == [2, 0]);

% 2D, create zonotope
c = [1; 5]; G = [2 3 4; 6 7 8];
Z = zonotope(c,G);
res(end+1,1) = all(withinTol(c,center(Z)));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
