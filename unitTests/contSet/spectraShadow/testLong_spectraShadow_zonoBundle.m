function res = testLong_spectraShadow_zonoBundle
% testLong_spectraShadow_zonoBundle - unit test function of zonoBundle
%
% Syntax:
%    res = testLong_spectraShadow_zonoBundle
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


% empty spectrahedron
SpS_empty = spectraShadow.empty(0);
zB_empty = zonoBundle(SpS_empty);
assert(isemptyobject(zB_empty));


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
zB = zonoBundle(SpS);
zB_true = zonoBundle(interval(-1,1));
assert(isequal(zB,zB_true,1e-10));

% 1D, empty
SpS = spectraShadow([-1 0]);
zB = zonoBundle(SpS);
zB_true = zonoBundle.empty(1);
assert(isequal(zB,zB_true,1e-10));

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
zB = zonoBundle(SpS);
zB_true = zonoBundle(interval(1,1));
assert(isequal(zB,zB_true,1e-10));


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
zB = zonoBundle(SpS);
zB_true = zonoBundle(interval([-1;-1],[1;1]));
assert(isequal(zB,zB_true,1e-10));

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
zB = zonoBundle(SpS);
zB_true = zonoBundle(interval([1;1],[1;1]));
assert(isequal(zB,zB_true,1e-10));

% 2D, empty
SpS = spectraShadow([-1 0 0]);
zB = zonoBundle(SpS);
zB_true = zonoBundle.empty(2);
assert(isequal(zB,zB_true,1e-10));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
