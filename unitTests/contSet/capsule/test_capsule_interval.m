function res = test_capsule_interval
% test_capsule_interval - unit test function of interval conversion
%
% Syntax:
%    res = test_capsule_interval
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

res = true(0);

% empty capsule
n = 2;
C = capsule.empty(n);
I = interval(C);
res(end+1,1) = representsa_(I,'emptySet',eps) && dim(I) == n;

% full-dimensional capsule
c = [2; 0; -1]; g = [0.2; -0.7; 0.4]; r = 1;
C = capsule(c,g,r);
% convert to interval
I = interval(C);
% check containment
res(end+1,1) = contains(I,C);

% one-dimensional capsule
c = 2; g = 1; r = 0;
C = capsule(c,g,r);
I = interval(C);
res(end+1,1) = isequal(I,interval(1,3));

% two-dimensional capsule with axis-aligned generator and no radius
c = [1;-1]; g = [1;0]; r = 0;
C = capsule(c,g,r);
I = interval(C);
res(end+1,1) = isequal(I,interval([0;-1],[2;-1]));

% two-dimensional capsule with all-zero generator and (no) radius
c = [0;-1]; g = [0;0]; r = 1;
C = capsule(c,g,r);
I = interval(C);
res(end+1,1) = isequal(I,interval([-1;-2],[1;0]));
c = [0;-1]; g = [0;0]; r = 0;
C = capsule(c,g,r);
I = interval(C);
res(end+1,1) = isequal(I,interval([0;-1]));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
