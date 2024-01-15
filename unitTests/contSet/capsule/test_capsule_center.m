function res = test_capsule_center
% test_capsule_center - unit test function of center
%
% Syntax:
%    res = test_capsule_center
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D empty capsule
n = 2;
C = capsule.empty(n);
c = center(C);
res(end+1,1) = isempty(c) && all(size(c) == [2, 0]);
    
% 3D capsule
c_true = [2; 0; -1]; g = [1; -1; 2]; r = 0.5;
C = capsule(c_true,g,r);
c = center(C);
res(end+1,1) = all(withinTol(c,c_true));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
