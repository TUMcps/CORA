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
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger, Tobias Ladner
% Written:       25-April-2022
% Last update:   27-July-2023 (MW, more cases)
%                31-October-2024 (TL, v-polytope)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% H-polytope --------------------------------------------------------------

% 1D, only inequalities, bounded
A = [2;-1]; b = [6;1];
P = polytope(A,b);
c = center(P);
c_true = 1;
assert(all(withinTol(c,c_true)));

% 1D, only equalities, single point
Ae = 3; be = 5;
P = polytope([],[],Ae,be);
c = center(P);
c_true = 5/3;
assert(all(withinTol(c,c_true)));

% 1D, only inequalities, unbounded
A = [3;2;4]; b = [5;2;-3];
P = polytope(A,b);
c = center(P);
assert(all(isnan(c)));

% 1D, only inequalities, empty
Ae = [1;4]; be = [2;-5];
P = polytope([],[],Ae,be);
c = center(P);
assert(isempty(c));

% 1D, inequalities and equalities, empty
A = [1;-4]; b = [4;-2]; Ae = 5; be = 100;
P = polytope(A,b,Ae,be);
c = center(P);
assert(isempty(c));

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
c = center(P);
assert(isscalar(c) && isnan(c));


% 2D, only inequalities, bounded
A = [1 1; -1 1; 1 -1; -1 -1]; b = ones(4,1);
P = polytope(A,b);
c = center(P);
c_true = [0; 0];
assert(all(withinTol(c,c_true)));

% 2D, only inequalities, empty
A = [1 0; -1 0]; b = [-1; -1];
P = polytope(A,b);
c = center(P);
assert(isempty(c));

% 2D, only equalities, empty
Ae = [1 0; 0 1; 0 1]; be = [1 -1 0];
P = polytope([],[],Ae,be);
c = center(P);
assert(isempty(c));

% 2D, only equalities, single point
Ae = [1 0; 0 1]; be = [0;0];
P = polytope([],[],Ae,be);
c = center(P);
c_true = [0; 0];
assert(all(withinTol(c,c_true)));

% 2D, inequalities and equalities, unbounded
A = [1 0]; b = 1; Ae = [0 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
assert(all(isnan(c)));

% 2D, inequalities and equalities, bounded
A = [1 0; -1 0]; b = [1; 1]; Ae = [0 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
c_true = [0; 1];
assert(all(withinTol(c,c_true)));

% 2D, inequalities and equalities, empty
A = [2 1; -1 2; 0 -1]; b = ones(3,1); Ae = [1 1]; be = 10;
P = polytope(A,b,Ae,be);
c = center(P);
assert(isempty(c));

% 2D, fully empty
Ae = zeros(0,2); be = zeros(0,0);
P = polytope([],[],Ae,be);
c = center(P);
assert(all(size(c) == [2,1]) && all(isnan(c)));

% 2D, V-polytope
V = [1 1; -1 1; -1 -1; 1 -1]';
P = polytope(V);
c = center(P);
c_true = [0;0];
assert(all(withinTol(c,c_true)));


% 3D, only inequalities, bounded
A = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1; 1 0 0; -1 0 0]; b = ones(6,1);
P = polytope(A,b);
c = center(P);
c_true = [0; 0; 0];
assert(all(withinTol(c,c_true)));

% 3D, inequalities and equalities, bounded, degenerate
A = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1]; b = ones(4,1); Ae = [1 0 0]; be = 2;
P = polytope(A,b,Ae,be);
c = center(P);
c_true = [2; 0; 0];
assert(all(withinTol(c,c_true)));

% 3D, only inequalities, unbounded
A = [1 0 0; 0 1 0]; b = [0; 0];
P = polytope(A,b);
c = center(P);
assert(all(isnan(c)));

% 3D, only equalities, unbounded
Ae = [1 0 0; 0 1 0]; be = [0;0];
P = polytope([],[],Ae,be);
c = center(P);
assert(all(isnan(c)));

% 3D, inequalities and equalities, unbounded
A = [1 1 0; -1 0 0]; b = [1; 1]; Ae = [0 1 1]; be = 1;
P = polytope(A,b,Ae,be);
c = center(P);
assert(all(isnan(c)));

% 3D, single point
Ae = [1 2 -1; 0 1 1; -1 2 1]; be = [1;1;1];
P = polytope([],[],Ae,be);
c = center(P);
c_true = [0.5;0.5;0.5];
assert(all(withinTol(c,c_true)));

% V-polytope --------------------------------------------------------------

% line
V = [1 3;2 4];
P = polytope(V);
c =center(P,'avg');
c_true = [2;3];
assert(all(withinTol(c,c_true),'all'));

% 1d
V = [1,3];
P = polytope(V);
c =center(P,'avg');
c_true = [2];
assert(isequal(c,c_true));

% 2d
V = [1 0 1; 0 1 1];
P = polytope(V);
c =center(P,'avg');
c_true = [2/3;2/3];
assert(all(withinTol(c,c_true),"all"));

% 2d
V = [1 0 1; 0 1 1];
P = polytope(V);
c =center(P,'chebyshev');
c_true = [ 0.7071067811865472 ; 0.7071067811865475 ];
assert(all(withinTol(c,c_true),"all"));

% 3d 
V = [1 0 1 1; 0 1 1 1; 0 0 0 1];
P = polytope(V);
c =center(P,'avg');
c_true = [0.75;0.75;0.25];
assert(all(withinTol(c,c_true),"all"));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
