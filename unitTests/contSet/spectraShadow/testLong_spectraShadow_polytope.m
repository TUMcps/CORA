function res = testLong_spectraShadow_polytope
% testLong_spectraShadow_polytope - unit test function of contains
%
% Syntax:
%    res = testLong_spectraShadow_polytope
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2 dimensional box with radius 3 around point [-1;2]:
A0 = eye(3);
A1 = [0 1 0;1 0 0;0 0 0];
A2 = [0 0 1;0 0 0;1 0 0];
SpS = spectraShadow([A0,A1,A2],[-1;2],3*eye(2));

% convert to polytope (inner-approximation)
P = polytope(SpS,'inner');
xs = [P.randPoint(100),P.randPoint(100,'extreme')];
assert(all(contains(SpS,xs,'exact',1e-8)));

% convert to polytope (outer-approximation)
P = polytope(SpS,'outer');
xs = [SpS.randPoint(100),SpS.randPoint(100,'extreme')];
assert(all(contains(P,xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
