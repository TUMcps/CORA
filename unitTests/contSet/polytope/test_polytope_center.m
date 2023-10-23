function res = test_polytope_center
% test_polytope_center - unit test function of center
%
% Syntax:
%    res = test_polytope_center
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

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       25-April-2022
% Last update:   27-July-2023 (MW, more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, only inequalities, bounded
P = polytope([2;-1],[6;1]);
c = center(P);
c_true = 1;
res(end+1,1) = all(withinTol(c,c_true));

% 1D, only equalities, single point
P = polytope([],[],3,5);
c = center(P);
c_true = 5/3;
res(end+1,1) = all(withinTol(c,c_true));

% 1D, only inequalities, unbounded
P = polytope([3;2;4],[5;2;-3]);
c = center(P);
res(end+1,1) = all(isnan(c));

% 1D, only inequalities, empty
P = polytope([],[],[1;4],[2;-5]);
c = center(P);
res(end+1,1) = isempty(c);

% 1D, inequalities and equalities, empty
P = polytope([1;-4],[4;-2],5,100);
c = center(P);
res(end+1,1) = isempty(c);


% 2D, only inequalities, bounded
P = polytope([1 1; -1 1; 1 -1; -1 -1],ones(4,1));
c = center(P);
c_true = [0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, only inequalities, empty
P = polytope([1 0; -1 0],[-1; -1]);
c = center(P);
res(end+1,1) = isempty(c);

% 2D, only equalities, empty
P = polytope([],[],[1 0; 0 1; 0 1],[1 -1 0]);
c = center(P);
res(end+1,1) = isempty(c);

% 2D, only equalities, single point
P = polytope([],[],[1 0; 0 1],[0;0]);
c = center(P);
c_true = [0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, inequalities and equalities, unbounded
P = polytope([1 0],1,[0 1],1);
c = center(P);
res(end+1,1) = all(isnan(c));

% 2D, inequalities and equalities, bounded
P = polytope([1 0; -1 0],[1; 1],[0 1],1);
c = center(P);
c_true = [0; 1];
res(end+1,1) = all(withinTol(c,c_true));

% 2D, inequalities and equalities, empty
P = polytope([2 1; -1 2; 0 -1],ones(3,1),[1 1],10);
c = center(P);
res(end+1,1) = isempty(c);


% 3D, only inequalities, bounded
P = polytope([0 1 0; 0 0 1; 0 -1 0; 0 0 -1; 1 0 0; -1 0 0],ones(6,1));
c = center(P);
c_true = [0; 0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 3D, inequalities and equalities, bounded, degenerate
P = polytope([0 1 0; 0 0 1; 0 -1 0; 0 0 -1],ones(4,1),[1 0 0],2);
c = center(P);
c_true = [2; 0; 0];
res(end+1,1) = all(withinTol(c,c_true));

% 3D, only inequalities, unbounded
P = polytope([1 0 0; 0 1 0],[0; 0]);
c = center(P);
res(end+1,1) = all(isnan(c));

% 3D, only equalities, unbounded
P = polytope([],[],[1 0 0; 0 1 0],[0;0]);
c = center(P);
res(end+1,1) = all(isnan(c));

% 3D, inequalities and equalities, unbounded
P = polytope([1 1 0; -1 0 0],[1; 1],[0 1 1],1);
c = center(P);
res(end+1,1) = all(isnan(c));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
