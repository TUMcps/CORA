function res = testLong_spectraShadow_polyZonotope
% testLong_spectraShadow_polyZonotope - unit test function of polyZonotope
%
% Syntax:
%    res = testLong_spectraShadow_polyZonotope
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

% We don't check for the exact execution of the algorithms, since there is
% no exact algorithm to check if two polynomial zonotopes are equal

% empty spectrahedron
SpS_empty = spectraShadow.empty();
pZ_empty = polyZonotope(SpS_empty);
assert(isemptyobject(pZ_empty));


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
pZ = polyZonotope(SpS);

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
pZ = polyZonotope(SpS);


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
pZ = polyZonotope(SpS);

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
pZ = polyZonotope(SpS);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
