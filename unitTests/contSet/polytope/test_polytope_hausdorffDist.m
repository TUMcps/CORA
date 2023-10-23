function res = test_polytope_hausdorffDist
% test_polytope_hausdorffDist - unit test function of hausdorffDist
%
% Syntax:
%    res = test_polytope_hausdorffDist
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
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% init polytope
A = [1 0; 0 1; -1 0; 0 -1];
b = ones(4,1);
P = polytope(A,b);

% point
p = [2;0];

% compute Hausdorff distance
dH = hausdorffDist(P,p);

% check result
dH_true = 1;
res(end+1,1) = withinTol(dH,dH_true);


% different point
p = [2; 2];
dH = hausdorffDist(P,p);

% check result
dH_true = sqrt(2);
res(end+1,1) = withinTol(dH,dH_true);


% point inside the polytope
p = [0;0];
dH = hausdorffDist(P,p);

% check result
dH_true = 0;
res(end+1,1) = withinTol(dH,dH_true);


% more complicated shape
A = [1 1; -1 1; -1 -1; 1 -1];
b = ones(4,1);
P = polytope(A,b);

p = [1;1];
dH = hausdorffDist(P,p);

% check result
dH_true = sqrt(2)/2;
res(end+1,1) = withinTol(dH,dH_true);


% multiple points (same result)
p = [1 1; 1 0; -0.5 1; 0 -0.5]';
dH = hausdorffDist(P,p);

% check result
dH_true = sqrt(2)/2;
res(end+1,1) = withinTol(dH,dH_true);

% degenerate case
D1 = polytope([1 0; 0 1;-1 0;0 -1],[1;1;1;-1]);
D2 = polytope([1 0;0 1;-1 0;0 -1],[6;1 ;-4; -1]);
dH = hausdorffDist(D1,D2);

% check result
dH_true = 5;
res(end+1,1) = withinTol(dH,dH_true);


% unbounded case
U = polytope([1 0; 0 1;-1 0],[1;1;1]);
p = [2;-100];
dH = hausdorffDist(U,p);

% check result
dH_true = 1.0;
res(end+1,1) = withinTol(dH,dH_true,1e-4);


% combine results
res
res = all(res);

% ------------------------------ END OF CODE ------------------------------
