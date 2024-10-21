function res = test_interval_radius
% test_interval_radius - unit test function of radius
%
% Syntax:
%    res = test_interval_radius
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
% Written:       27-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% create interval (all positive)
Int = interval([-2;-3], [4;1]);

% compute radius
Int_rad = radius(Int);

% true solution
true_rad = sqrt(13);
assert(Int_rad == true_rad);

% matrix
I = interval([-1.5;1.1],[1.2;2.4]);
true_rad = 1.498332406377170;
assert(withinTol(true_rad,radius(I)))

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
true_rad = 3.889087296526011;
assert(withinTol(radius(I),true_rad))

% ------------------------------ END OF CODE ------------------------------
