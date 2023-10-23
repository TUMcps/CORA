function res = test_polytope_volume
% test_polytope_volume - unit test function of volume
%
% Syntax:
%    res = test_polytope_volume
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
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% empty case
P = polytope();

% compute volume
val = volume(P);

% true volume
val_true = 0;
res(end+1,1) = withinTol(val,val_true);

% 1D polytopes (bounded and unbounded)
P = polytope(1,0);
val = volume(P);
res(end+1,1) = val == Inf;
P = polytope([1;-2],[2;6]);
val = volume(P);
res(end+1,1) = withinTol(val,5);


% init polytope
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);

% compute volume
val = volume(P);

% true volume (computed by hand)
val_true = 2*2 + 1*2 + 0.5*2*2 + 0.5*3*3 + 2*2 + 0.5*1*3 + 0.5*1*4;
res(end+1,1) = withinTol(val,val_true);


% volume should be unaffected by translation
z = [2; 1];
V = V - z;
P = polytope(V);

% compute volume
val = volume(P);
res(end+1,1) = withinTol(val,val_true);


% degenerate polytope
A = [1 1 0; -1 1 0; 0 -1 0; 0 0 1; 0 0 -1];
b = [1; 1; 1; 0; 0];
P = polytope(A,b);

% compute volume
val = volume(P);
res(end+1,1) = withinTol(val,0);


% unbounded polytope
A = [1 1 0; -1 1 0; 0 -1 0];
b = ones(3,1);
P = polytope(A,b);

% compute volume
val = volume(P);
res(end+1,1) = val == Inf;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
