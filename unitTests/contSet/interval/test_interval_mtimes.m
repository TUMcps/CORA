function res = test_interval_mtimes
% test_interval_mtimes - unit test function of mtimes,
%    overloaded '*' operator for intervals
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
% See also: -

% Authors:       Matthias Althoff
% Written:       05-August-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%init result vector
resVec = [];

% empty set
a = interval.empty(1);
b = 3;
resVec(end+1) = representsa(b*a,'emptySet');

%SCALAR VALUES-------------------------------------------------------------
%test 1: a: interval, b: numeric
a = interval(-1,2);
b = 3;
c = a*b;
c_true = interval(-3,6);
resVec(end+1) = (c == c_true);

%test 2: a: numeric, b: interval
a = -2;
b = interval(2,4);
c = a*b;
c_true = interval(-8,-4);
resVec(end+1) = (c == c_true);

%test 3: a: interval, b: interval
a = interval(-2,1);
b = interval(2,4);
c = a*b;
c_true = interval(-8,4);
resVec(end+1) = (c == c_true);
%--------------------------------------------------------------------------

%MIXED SCALAR/MATRIX VALUES------------------------------------------------
%a: scalar, b: matrix 
%test 4: a: interval, b: numeric
a = interval(-1,2);
b = [-1 0; 1 2];
c = a*b;
c_true = interval([-2 0; -1 -2],[1 0; 2 4]);
resVec(end+1) = (c == c_true);

%test 5: a: numeric, b: interval
a = -2;
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-8 4; -2 -6],[-4 6; 2 -4]);
resVec(end+1) = (c == c_true);

%test 6: a: interval, b: interval
a = interval(-1,2);
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-4 -6; -2 -3],[8 3; 2 6]);
resVec(end+1) = (c == c_true);

%a: matrix, b: scalar
%test 7: a: interval, b: numeric
a = interval([-1 0; -2 2],[2 1; -1 3]);
b = -1;
c = a*b;
c_true = interval([-2 -1; 1 -3],[1 0; 2 -2]);
resVec(end+1) = (c == c_true);

%test 8: a: numeric, b: interval
a = [-1 0; 1 2];
b = interval(-2,1);
c = a*b;
c_true = interval([-1 0; -2 -4],[2 0; 1 2]);
resVec(end+1) = (c == c_true);

%test 9: a: interval, b: interval
a = interval([-1 0; -2 2],[2 1; -1 3]);
b = interval(-2,1);
c = a*b;
c_true = interval([-4 -2; -2 -6],[2 1; 4 3]);
resVec(end+1) = (c == c_true);
%--------------------------------------------------------------------------

%MATRIX VALUES-------------------------------------------------------------
%test 10: a: interval, b: numeric
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = [-1 0; 1 2];
c = a*b;
c_true = interval([-7 -6; 1 4],[-4 -4; 4 6]);
resVec(end+1) = (c == c_true);

%test 11: a: numeric, b: interval
a = [-2 1; -3 2];
b = interval([2 -3; -1 2],[4 -2; 1 3]);
c = a*b;
c_true = interval([-9 6; -14 10],[-3 9; -4 15]);
resVec(end+1) = (c == c_true);

%test 12: a: interval, b: interval
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = interval([-2 0; -1 2],[-1 1; 1 3]);
c = a*b;
c_true = interval([-11 -9; -5 3],[1 0; 5 10]);
resVec(end+1) = (c == c_true);
%--------------------------------------------------------------------------

%SPARSE VALUES-------------------------------------------------------------
a = interval([2 -3; -1 2],[4 -2; 1 3]);
b = interval(sparse([-2 0; -1 2]),sparse([-1 1; 1 3]));
c = a*b;
c_true = interval([-11 -9; -5 3],[1 0; 5 10]);
resVec(end+1) = (c == c_true);

% Unbounded intervals and 0 -----------------------------------------------

% scalar
I0 = interval(0);
I = interval(-inf,inf);
resVec(end+1) = isequal(I0,0*I);

% higher-dim
I0 = interval(zeros(3,1));
I = interval(-inf(3,1),inf(3,1));
resVec(end+1) = isequal(I0,0*I);

% -------------------------------------------------------------------------

% check result
res = all(resVec);

% ------------------------------ END OF CODE ------------------------------
