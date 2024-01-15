function res = test_polytope_setProperty_minHRep
% test_polytope_setProperty_minHRep - unit test function to check whether
%    the internally-used set property 'minHRep' is changed correctly
%    following different set operations on a polytope
%
% Syntax:
%    res = test_polytope_setProperty_minHRep
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
% Written:       01-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% --- box -----------------------------------------------------------------

% 2D, empty set
P = polytope([1 0; -1 0],[2;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;

% 2D, unbounded, non-degenerate, non-empty
P = polytope([1 0 0; 0 1 0],[1;-3]);
P_ = box(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;

% 2D, bounded, degenerate, non-empty
P = polytope([1 0 0; 0 1 0; -1 -1 0],[1;2;2],[0,0,1],0);
P_ = box(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;


% --- compact -------------------------------------------------------------

% 1D, empty
P = polytope([1;-1],[1;-2]);
P_ = compact(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;
% 3D, empty
P = polytope([1 1 0; -1 1 0; 0 -1 0; 1 0 0],[1;1;1;-3],[0 0 1],3);
P_ = compact(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;


% --- empty ---------------------------------------------------------------

n = 2;
P = polytope.empty(n);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;


% --- Inf -----------------------------------------------------------------

n = 2;
P = polytope.Inf(n);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;


% --- plus ----------------------------------------------------------------

% 2D, bounded + vector
A = [1 0; -1 1; -1 -1]; b = ones(3,1);
P = polytope(A,b);
compact(P);
v = [-1;1];
P_sum = P + v;
% resulting polytope is also non-empty
res(end+1,1) = ~isempty(P_sum.minHRep.val) && P_sum.minHRep.val;


% --- polytope ------------------------------------------------------------

% 2D, only inequalities, non-empty
P = polytope([1 1; -1 1; 0 -1],ones(3,1));
% remove redundancies
P = compact(P);
% copy polytope, property should also be copied
P_ = polytope(P);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;
% 1D, with redundancies
P = polytope([2;1],[1;1]);
res(end+1,1) = isempty(P.minHRep.val);
% init set via vertex representation
% 1D, unbounded, non-degenerate
V = [-Inf, 2];
P = polytope(V);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
% 1D, bounded, non-degenerate
V = [-3 -1 4 5];
P = polytope(V);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
% 1D, bounded, degenerate
V = 2;
P = polytope(V);
res(end+1,1) = ~isempty(P.minHRep.val) && P.minHRep.val;
% 2D, bounded, non-degenerate
V = [2 1; -1 4; -4 0; -1 -2; 3 -1]';
P = polytope(V);
% redundancy unknown
res(end+1,1) = isempty(P.minHRep.val);
% 2D, bounded, degenerate
V = [-1 1; 2 0]';
P = polytope(V);
% redundancy unknown
res(end+1,1) = isempty(P.minHRep.val);


% --- project -------------------------------------------------------------

% 3D, empty
P = polytope([1 0 0; -1 0 0],[2;-3]);
P_ = project(P,[1,2]);
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;


% --- lift ----------------------------------------------------------------

% 2D, non-empty
P = polytope([1 1; -1 1; 0 -1],[1;1;1]);
% remove redundancies
P = compact(P);
% project to higher-dimensional space
P_ = lift(P,5,[2,3]);
% higher-dimensional polytope also non-empty
res(end+1,1) = ~isempty(P_.minHRep.val) && P_.minHRep.val;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
