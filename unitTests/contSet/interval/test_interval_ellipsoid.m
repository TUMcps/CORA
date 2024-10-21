function res = test_interval_ellipsoid
% test_interval_ellipsoid - unit test function of interval/ellipsoid
%
% Syntax:
%    res = test_interval_ellipsoid
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

% Authors:       Tobias Ladner
% Written:       18-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check outer approximation ---

I = interval([-1; -2; 1], [3; 1; 4]);
E = ellipsoid(I);
assert(contains(E,I));

E = ellipsoid(I, 'outer');
assert(contains(E,I));

% check inner approximation ---

E = ellipsoid(I, 'inner');
assert(contains(I,E));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
