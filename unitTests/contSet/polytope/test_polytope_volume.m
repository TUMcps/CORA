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

res = true(0);

% empty case
P = polytope.empty(2);
val = volume(P);
val_true = 0;
res(end+1,1) = withinTol(val,val_true);

% 1D, bounded
A = [1;-2]; b = [2;6];
P = polytope(A,b);
val = volume(P);
val_true = 5;
res(end+1,1) = withinTol(val,val_true);

% 1D, unbounded
A = 1; b = 0;
P = polytope(A,b);
val = volume(P);
val_true = Inf;
res(end+1,1) = val == val_true;


% 2D, fully empty (unbounded)
A = zeros(0,2); b = zeros(0,0);
P = polytope(A,b);
val = volume(P);
val_true = Inf;
res(end+1,1) = val == val_true;

% 2D, vertex instantiation
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
val = volume(P);
% true volume (computed by hand)
val_true = 2*2 + 1*2 + 0.5*2*2 + 0.5*3*3 + 2*2 + 0.5*1*3 + 0.5*1*4;
res(end+1,1) = withinTol(val,val_true);
% translate vertices (should not affect the volume)
V = V - [2; 1];
P = polytope(V);
val = volume(P);
res(end+1,1) = withinTol(val,val_true);

% 2D, unbounded
A = [1 0; 1 -1; -1 -1]; b = ones(3,1);
P = polytope(A,b);
val = volume(P);
res(end+1,1) = val == Inf;


% 3D, degenerate polytope
A = [1 1 0; -1 1 0; 0 -1 0; 0 0 1; 0 0 -1]; b = [1; 1; 1; 0; 0];
P = polytope(A,b);
val = volume(P);
res(end+1,1) = withinTol(val,0);

% 3D, unbounded polytope
A = [1 1 0; -1 1 0; 0 -1 0]; b = ones(3,1);
P = polytope(A,b);
val = volume(P);
res(end+1,1) = val == Inf;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
