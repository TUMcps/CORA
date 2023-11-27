function res = test_capsule_polytope
% test_capsule_polytope - unit test function for the conversion to polytope
%
% Syntax:
%    res = test_capsule_polytope
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
% Written:       25-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, just center
C = capsule(2);
P = polytope(C);
P_ = polytope([1;-1],[2;-2]);
res(end+1,1) = P == P_;

% 1D, center and generator
C = capsule(1,-1);
P = polytope(C);
P_ = polytope([1;-1],[2;0]);
res(end+1,1) = P == P_;

% 1D, 'full' capsule
C = capsule(1,-1,0.5);
P = polytope(C);
P_ = polytope([1;-1],[2.5;0.5]);
res(end+1,1) = P == P_;


% 2D, just center
C = capsule([1;-1]);
P = polytope(C);
P_ = polytope([1;-1]);
res(end+1,1) = P == P_;

% 2D, just a line
C = capsule([1;-1],[4;3],0);
P = polytope(C);
P_ = polytope([4/5 3/5; -4/5 -3/5; -3 4; 3 -4],[5;5;0;0]) + [1;-1];
res(end+1,1) = P == P_;

% 2D, 'full' capsule
C = capsule([0;0],[4;3],2);
P = polytope(C,'outer');
res(end+1,1) = contains(P,C);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
