function res = test_interval_display
% test_interval_display - unit test function of display
%
% Syntax:
%    res = test_interval_display
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
% Written:       21-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty interval
I = interval.empty(2)

% column interval
I = interval([-3;2;1],[4;5;6])

% row interval
I = interval([1 2 3],[4 5 6]);

% double values
I = interval([-0.253 -0.128 -0.654],[0.752 0.983 0.127])

% matrices
I = interval(-magic(3),magic(3))

% sparse
lb = [-1 0 -2; 0 0 -3; 0 0 0];
ub = [2 1 0; 0 0 0; 0 0 1];
I = interval(sparse(lb),ub)
I = interval(lb,sparse(ub))
I = interval(sparse(lb),sparse(ub))

% sparse with all-zero entries
lb = [0 0 0; 0 0 0; 0 0 0];
ub = [0 0 0; 0 0 0; 0 0 0];
I = interval(sparse(lb),ub)
I = interval(lb,sparse(ub))
I = interval(sparse(lb),sparse(ub))

% large sparse matrix
lb = [-1; zeros(14,1)];
ub = [zeros(14,1); 1];
I = interval(sparse(lb),sparse(ub))

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub)

% ------------------------------ END OF CODE ------------------------------
