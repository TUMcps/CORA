function res = test_spectraShadow_dim
% test_spectraShadow_dim - unit test function of dim
%
% Syntax:
%    res = test_spectraShadow_dim
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
assert(dim(SpS_empty) == 0);


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
assert(dim(SpS) == 1);

% 1D, empty
SpS = spectraShadow([-1 0]);
assert(dim(SpS) == 1);

% 1D, unbounded
SpS = spectraShadow([1 0]);
assert(dim(SpS) == 1);

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
assert(dim(SpS) == 1);


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(dim(SpS) == 2);

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(dim(SpS) == 2);

% 2D, unbounded, non-degenerate 
SpS = spectraShadow([1 0 0]);
assert(dim(SpS) == 2);

% 2D, unbounded, degenerate
A0 = [-1 0;0 1];
Ai{1} = zeros(2);
Ai{2} = [1 0;0 -1];
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
assert(dim(SpS) == 2);

% 2D, empty
SpS = spectraShadow([-1 0 0]);
assert(dim(SpS) == 2);


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
