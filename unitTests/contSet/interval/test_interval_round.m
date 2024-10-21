function res = test_interval_round
% test_interval_round - unit test function of radius
%
% Syntax:
%    res = test_interval_round
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

% Authors:       Tobias Ladner
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% scalar
I = interval(1.42,1.71);
I_rounded = round(I);
assert(isequal(I_rounded,interval(1,2)))
I_rounded = round(I,1);
assert(isequal(I_rounded,interval(1.4,1.7)))

% matrix
inf = [ -1.610 -0.766 ; -0.335 -2.862 ];
sup = [ -1.221 1.771 ; 1.886 -0.788 ];
I = interval(inf,sup);
inf = [ -1.600 -0.800 ; -0.300 -2.900 ];
sup = [ -1.200 1.800 ; 1.900 -0.800 ];
I_rounded = interval(inf,sup);
assert(isequal(round(I,1),I_rounded))

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
assert(isequal(round(I),interval(round(lb),round(ub))))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
