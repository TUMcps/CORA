function res = test_capsule_isFullDim
% test_capsule_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_capsule_isFullDim
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% define properties
c = [2; 0; -1];
g = [-1; 1; 2];
g0 = [0; 0; 0];
r = 0.5;
r0 = 0;

% generator and radius all-zero
C = capsule(c,g0,r0);
assert(~isFullDim(C));

% generator all-zero
C = capsule(c,g0,r);
assert(isFullDim(C));

% radius is zero
C = capsule(c,g,r0);
assert(~isFullDim(C));

% generator and radius non-zero
C = capsule(c,g,r);
assert(isFullDim(C));

% empty set
C = capsule.empty(2);
assert(~isFullDim(C));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
