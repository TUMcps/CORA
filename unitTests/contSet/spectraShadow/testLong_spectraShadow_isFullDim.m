function res = testLong_spectraShadow_isFullDim
% testLong_spectraShadow_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_spectraShadow_isFullDim
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
assert(~isFullDim(SpS_empty));


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
assert(isFullDim(SpS));

% 1D, empty
SpS = spectraShadow([-1 0]);
assert(~isFullDim(SpS));

% 1D, unbounded
SpS = spectraShadow([1 0]);
assert(isFullDim(SpS));

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
assert(~isFullDim(SpS));


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(isFullDim(SpS));

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(~isFullDim(SpS));

% 2D, unbounded, non-degenerate 
SpS = spectraShadow([1 0 0]);
assert(isFullDim(SpS));

% 2D, unbounded, degenerate
A0 = [-1 0;0 1];
Ai{1} = zeros(2);
Ai{2} = [1 0;0 -1];
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(~isFullDim(SpS));

% 2D, empty
SpS = spectraShadow([-1 0 0]);
assert(~isFullDim(SpS));
% Check if fullDim hidden property is set
assert(~isempty(SpS.fullDim.val) && SpS.fullDim.val == false);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
