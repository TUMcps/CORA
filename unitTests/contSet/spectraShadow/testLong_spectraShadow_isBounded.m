function res = testLong_spectraShadow_isBounded
% testLong_spectraShadow_isBounded - unit test function of isBounded
%
% Syntax:
%    res = testLong_spectraShadow_isBounded
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
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% empty spectrahedron
SpS_empty = spectraShadow.empty();
assert(isBounded(SpS_empty) == 1);


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
assert(isBounded(SpS) == 1);

% 1D, empty
SpS = spectraShadow([-1 0]);
assert(isBounded(SpS) == 1);

% 1D, unbounded
SpS = spectraShadow([1 0]);
assert(isBounded(SpS) == 0);

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
assert(isBounded(SpS) == 1);


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(isBounded(SpS) == 1);

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(isBounded(SpS) == 1);

% 2D, unbounded, non-degenerate 
SpS = spectraShadow([1 0 0]);
assert(isBounded(SpS) == 0);

% 2D, unbounded, degenerate
A0 = [-1 0;0 1];
Ai{1} = zeros(2);
Ai{2} = [1 0;0 -1];
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(isBounded(SpS) == 0);

% 2D, empty
SpS = spectraShadow([-1 0 0]);
assert(isBounded(SpS) == 1);

% Checking that the hidden property 'bounded' is properly set:
SpS = spectraShadow.generateRandom();
bounded = isBounded(SpS);
assert(bounded == SpS.bounded.val);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
