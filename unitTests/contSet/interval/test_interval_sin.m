function res = test_interval_sin
% test_interval_sin - unit_test_function of sine for intervals,
%    overloaded 'sin()' function for intervals
%
% Syntax:
%    res = test_interval_sin
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       13-January-2016
% Last update:   03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% bounded
I = interval([-0; 0.0; 0.0; 0; 0; -pi/4; -2*pi], ...
             [pi/4; pi/2; pi; 2*pi; 4*pi; pi/4; 4*pi]);
I_sin = sin(I);
I_true = interval([0;0;0;-1;-1;-sqrt(2)/2;-1],[sqrt(2)/2;1;1;1;1;sqrt(2)/2;1]);
assert(isequal(I_sin,I_true,tol));

% unbounded
I = interval([-Inf;0;-Inf],[Inf;Inf;0]);
I_sin = sin(I);
I_true = interval(-ones(3,1),ones(3,1));
assert(isequal(I_sin,I_true,tol));

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(inf,sup);
Isin = sin(I);
inf = reshape([ 0.841 -0.757 -0.757 -0.959 -1.000 0.000 0.141 0.841 0.000 -1.000 -0.841 -0.757 0.000 0.000 0.000 0.000 0.841 -0.841 0.141 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 0.997 0.141 0.909 1.000 -0.141 0.000 0.909 1.000 0.841 0.000 1.000 0.141 0.000 0.000 0.000 0.000 1.000 -0.479 1.000 1.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I_true = interval(inf,sup);
assert(isequal(Isin,I_true,1e-3))


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
