function res = test_interval_mtimes
% test_interval_mtimes - unit test function of overloaded '*' operator for
%    intervals, including uncertain linear maps of sets
%
% Syntax:
%    res = test_interval_mtimes
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

% Authors:       Matthias Althoff
% Written:       05-August-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
a = interval.empty(1);
b = 3;
assert(representsa(b*a,'emptySet'));

%SCALAR VALUES-------------------------------------------------------------
%test 1: a: interval, b: numeric
a = interval(-1,2);
b = 3;
c = a*b;
c_true = interval(-3,6);
assert(c == c_true);

%test 2: a: numeric, b: interval
a = -2;
b = interval(2,4);
c = a*b;
c_true = interval(-8,-4);
assert(c == c_true);

%test 3: a: interval, b: interval
a = interval(-2,1);
b = interval(2,4);
c = a*b;
c_true = interval(-8,4);
assert(c == c_true);
%--------------------------------------------------------------------------

%MIXED SCALAR/MATRIX VALUES------------------------------------------------
%a: scalar, b: matrix 
%test 4: a: interval, b: numeric
a = interval(-1,2);
b = [-1 0; 1 2];
c = a*b;
c_true = interval([-2 0; -1 -2],[1 0; 2 4]);
assert(c == c_true);

%test 5: a: numeric, b: interval
a = -2;
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-8 4; -2 -6],[-4 6; 2 -4]);
assert(c == c_true);

%test 6: a: interval, b: interval
a = interval(-1,2);
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-4 -6; -2 -3],[8 3; 2 6]);
assert(c == c_true);

%a: matrix, b: scalar
%test 7: a: interval, b: numeric
a = interval([-1 0; -2 2],[2 1; -1 3]);
b = -1;
c = a*b;
c_true = interval([-2 -1; 1 -3],[1 0; 2 -2]);
assert(c == c_true);

%test 8: a: numeric, b: interval
a = [-1 0; 1 2];
b = interval(-2,1);
c = a*b;
c_true = interval([-1 0; -2 -4],[2 0; 1 2]);
assert(c == c_true);

%test 9: a: interval, b: interval
a = interval([-1 0; -2 2],[2 1; -1 3]);
b = interval(-2,1);
c = a*b;
c_true = interval([-4 -2; -2 -6],[2 1; 4 3]);
assert(c == c_true);
%--------------------------------------------------------------------------

%MATRIX VALUES-------------------------------------------------------------
%test 10: a: interval, b: numeric
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = [-1 0; 1 2];
c = a*b;
c_true = interval([-7 -6; 1 4],[-4 -4; 4 6]);
assert(c == c_true);

%test 11: a: numeric, b: interval
a = [-2 1; -3 2];
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-9 6; -14 10],[-3 9; -4 15]);
assert(c == c_true);

%test 12: a: interval, b: interval
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = interval([-2 0; -1 2],[-1 1; 1 3]);
c = a*b;
c_true = interval([-11 -9; -5 3],[1 0; 5 10]);
assert(c == c_true);
%--------------------------------------------------------------------------

%SPARSE VALUES-------------------------------------------------------------
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = interval(sparse([-2 0; -1 2]),sparse([-1 1; 1 3]));
c = a*b;
c_true = interval([-11 -9; -5 3],[1 0; 5 10]);
assert(c == c_true);

% Unbounded intervals and 0 -----------------------------------------------

% scalar
I0 = interval(0);
I = interval(-inf,inf);
assert(isequal(I0,0*I));

% higher-dim
I0 = interval(zeros(3,1));
I = interval(-inf(3,1),inf(3,1));
assert(isequal(I0,0*I));

I = interval([1,2;-Inf,-Inf],[Inf,3;Inf,Inf]);
assert(isequal([0 0] * I, interval([0,0])))

% -------------------------------------------------------------------------

% intervals and contSet ---------------------------------------------------

% scalar interval
I = interval(0.5,2);

% ...zonotope
Z = zonotope([2;1],[1 0 -1; 0 1 1]);
Z_mtimes = I*Z;
assert(contains(Z_mtimes,I.inf*Z));
assert(contains(Z_mtimes,I.sup*Z));

% ...polyZonotope
c = [1; 2]; G = [1 2 1 -3; 1 -1 2 -1]; E = [1 0 0 2; 0 1 2 1];
pZ = polyZonotope(c,G,[],E);
pZ_mtimes = I*pZ;
c_true = [1.25; 2.5];
G_true =  [1.25 2.5 1.25 -3.75; 1.25 -1.25 2.5 -1.25];
GI_true = [6 0; 0 5.25];
assert(compareMatrices(pZ_mtimes.c,c_true));
assert(compareMatrices(pZ_mtimes.G,G_true));
assert(compareMatrices(pZ_mtimes.GI,GI_true));
assert(compareMatrices(pZ_mtimes.E,pZ.E));
assert(compareMatrices(pZ_mtimes.id,pZ.id));

% interval matrix
I = interval([1 -0.5; -1 0], [3 0.5; 3 2]);

% ...zonotope
c = [1; 2]; G = [1 2 1 -3; 1 -1 2 -1];
Z = zonotope(c,G);
Z_mtimes = I * Z;
assert(contains(Z_mtimes,I.inf*Z));
assert(contains(Z_mtimes,I.sup*Z));

% ...polyZonotope
c = [1; 2]; G = [1 2 1 -3; 1 -1 2 -1]; E = [1 0 0 2; 0 1 2 1];
pZ = polyZonotope(c,G,[],E);
pZ_mtimes = I * pZ;
c_true = [2; 3]; G_true =  [2 4 2 -6; 2 1 3 -4]; GI_true = [11.5 0; 0 23];
assert(compareMatrices(pZ_mtimes.c,c_true));
assert(compareMatrices(pZ_mtimes.G,G_true));
assert(compareMatrices(pZ_mtimes.GI,GI_true));
assert(compareMatrices(pZ_mtimes.E,pZ.E));
assert(compareMatrices(pZ_mtimes.id,pZ.id));

% n-d arrays (only w/ scalar)
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
I_mtimes = 2*I;
I_true = interval(2*lb,2*ub);
assert(isequal(I_mtimes,I_true));


% -------------------------------------------------------------------------

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
