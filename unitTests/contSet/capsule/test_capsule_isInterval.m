function res = test_capsule_isInterval
% test_capsule_isInterval - unit test function of isInterval
%
% Syntax:  
%    res = test_capsule_isInterval
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

% Author:       Mark Wetzlinger
% Written:      24-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty capsule
C = capsule();
res = isInterval(C);

% full-dimensional capsule
c = [2; 0; -1];
g = [0.2; -0.7; 0.4];
r = 1;
C = capsule(c,g,r);
res(end+1,1) = ~isInterval(C);

% one-dimensional capsule
C = capsule(2,1,0);
res(end+1,1) = isInterval(C);

% two-dimensional capsule with axis-aligned generator and no radius
C = capsule([1;-1],[1;0],0);
res(end+1,1) = isInterval(C);

% two-dimensional capsule with all-zero generator and (no) radius
C = capsule([0;-1],[0;0],1);
res(end+1,1) = ~isInterval(C);
C = capsule([0;-1],[0;0],0);
res(end+1,1) = isInterval(C);

% combine results
res = all(res);

%------------- END OF CODE --------------