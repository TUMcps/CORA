function res = test_spectraShadow_copy
% test_spectraShadow_copy - unit test function of copy
%
% Syntax:
%    res = test_spectraShadow_copy
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D spectraShadow
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
SpS_copy = copy(SpS);

% ensure that set properties are independent from another
isFullDim(SpS);
assert(isempty(SpS_copy.fullDim.val));

% ensure that set properties are copied
SpS_copy = copy(SpS);
assert(SpS_copy.fullDim.val);

% ensure that sets are equal (currently no isequal implemented)
% assert(isequal(SpS,SpS_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
