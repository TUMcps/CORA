function res = test_capsule_isnan
% test_capsule_isnan - unit test function of isnan
%
% Syntax:
%    res = test_capsule_isnan
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
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty capsule
C = capsule.empty(2);
assert(~isnan(C));

% 3D capsule
C = capsule([1;0;-1],[1;0;0],1);
assert(~isnan(C));

% no radius
C = capsule([1;0;-1],[1;0;0],0);
assert(~isnan(C));

% all-zero generator
C = capsule([1;0;-1],[0;0;0],0);
assert(~isnan(C));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
