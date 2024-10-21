function res = test_interval_enlarge
% test_interval_enlarge - unit test function of enlarge
%
% Syntax:
%    res = test_interval_enlarge
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
% Written:       29-August-2019
% Last update:   03-December-2023 (MW, add unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% bounded
I = interval([-2; -4; -3],[2; 3; 1]);
% ...bounded scaling factor
factor = 2;
I_enlarge = enlarge(I, factor);
I_true = interval([-4; -7.5; -5],[4; 6.5; 3]);
assert(isequal(I_enlarge,I_true));

% ...Inf as scaling factor
factor = Inf;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(3,1),Inf(3,1));
assert(isequal(I_enlarge,I_true));

% ...0 as scaling factor
factor = 0;
I_enlarge = enlarge(I, factor);
I_true = interval(center(I));
assert(isequal(I_enlarge,I_true));


% unbounded
I = interval([-Inf;-2],[2;Inf]);
% ...bounded scaling factor
factor = 2;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(2,1),Inf(2,1));
assert(isequal(I_enlarge,I_true));

% ...Inf as scaling factor
factor = Inf;
I_enlarge = enlarge(I, factor);
I_true = interval(-Inf(2,1),Inf(2,1));
assert(isequal(I_enlarge,I_true));

% ...0 as scaling factor (unbounded case throws an error)
factor = 0;
assertThrowsAs(@enlarge,'CORA:notSupported',I,factor);

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
I_enlarge = enlarge(I,2);
inf = reshape([ 0.750 2.500 1.000 2.500 -4.000 0.000 1.500 0.500 -0.500 -3.000 -2.500 2.500 0.000 0.000 0.000 0.000 0.500 -1.250 0.000 -1.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.750 4.500 5.000 12.500 0.000 0.000 3.500 2.500 1.500 1.000 3.500 4.500 0.000 0.000 0.000 0.000 2.500 -0.250 4.000 3.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I_true = interval(inf,sup);
assert(isequal(I_enlarge,I_true))

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
