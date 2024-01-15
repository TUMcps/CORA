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

res = true(0);

% define properties
c = [2; 0; -1];
g = [-1; 1; 2];
g0 = [0; 0; 0];
r = 0.5;
r0 = 0;

% generator and radius all-zero
C = capsule(c,g0,r0);
res(end+1,1) = ~isFullDim(C);

% generator all-zero
C = capsule(c,g0,r);
res(end+1,1) = isFullDim(C);

% radius is zero
C = capsule(c,g,r0);
res(end+1,1) = ~isFullDim(C);

% generator and radius non-zero
C = capsule(c,g,r);
res(end+1,1) = isFullDim(C);

% empty set
C = capsule.empty(2);
res(end+1,1) = ~isFullDim(C);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
