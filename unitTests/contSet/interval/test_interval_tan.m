function res = test_interval_tan
% test_interval_tan - unit test function of tangent function
%
% Syntax:
%    res = test_interval_tan
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
% See also: mtimes

% Authors:       Daniel Althoff, Mark Wetzlinger
% Written:       03-November-2015
% Last update:   04-January-2016
%                04-December-2023 (MW, update syntax, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% bounded, degenerate
I = interval(0);
I_tan = tan(I);
assert(isequal(I_tan,I,tol));

I = interval(1);
I_tan = tan(I);
I_true = interval(1.55740772465491);
assert(isequal(I_tan,I_true,tol));

% bounded, non-degenerate
I = interval(-1,0);
I_tan = tan(I);
I_true = interval(-1.55740772465491,0);
assert(isequal(I_tan,I_true,tol));

I = interval(-2,0);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
assert(isequal(I_tan,I_true,tol));

I = interval(-pi/2,0);
I_tan = tan(I);
I_true = interval(-1.63312393531954e+16,0);
assert(isequal(I_tan,I_true,tol));

I = interval(-pi/2,pi/2);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
assert(isequal(I_tan,I_true,tol));

% unbounded, non-degenerate
I = interval(-Inf,0);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
assert(isequal(I_tan,I_true,tol));

I = interval(-Inf,Inf);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
assert(isequal(I_tan,I_true,tol));

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(inf,sup);
Itan = tan(I);
inf = reshape([ 1.557 -0.143 -2.185 -Inf -Inf 0.000 -2.185 -Inf 0.000 -Inf -Inf -0.143 0.000 0.000 0.000 0.000 -Inf -1.557 -Inf -Inf 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 14.101 1.158 1.158  Inf  Inf 0.000 -0.143  Inf 1.557  Inf  Inf 1.158 0.000 0.000 0.000 0.000  Inf -0.546  Inf  Inf 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I_true = interval(inf,sup);
assert(isequal(Itan,I_true,1e-3))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
