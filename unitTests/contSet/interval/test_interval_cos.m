function res = test_interval_cos
% test_interval_cos - unit test function of cosine for intervals,
%    overloaded 'cos()' function for intervals
%
% Syntax:
%    res = test_interval_cos
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
% Written:       05-January-2016
% Last update:   03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% bounded
I = interval([-0; 0.0; 0.0; 0; 0; -pi/4; -2*pi], ...
             [pi/4; pi/2; pi; 2*pi; 4*pi; pi/4; 4*pi]);
I_cos = cos(I);
I_true = interval([sqrt(2)/2;0;-1;-1;-1;sqrt(2)/2;-1],ones(7,1));
assert(isequal(I_cos,I_true,tol));

% unbounded
I = interval([-Inf;0;-Inf],[Inf;Inf;0]);
I_cos = cos(I);
I_true = interval(-ones(3,1),ones(3,1));
assert(isequal(I_cos,I_true,tol));

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(inf,sup);
Icos = cos(I);
inf = reshape([ 0.071 -1.000 -1.000 -1.000 -0.990 1.000 -0.990 -0.416 0.540 -0.416 -0.416 -1.000 1.000 1.000 1.000 1.000 -0.416 0.540 -0.990 -0.416 1.000 1.000 1.000 1.000 ], [2,2,2,3]);
sup = reshape([ 0.540 -0.654 -0.416 1.000 0.540 1.000 -0.416 0.540 1.000 1.000 1.000 -0.654 1.000 1.000 1.000 1.000 0.540 0.878 0.540 1.000 1.000 1.000 1.000 1.000 ], [2,2,2,3]);
I_true = interval(inf,sup);
assert(isequal(Icos,I_true,1e-3))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
