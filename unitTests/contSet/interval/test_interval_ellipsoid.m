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

resvec = true(0);

% check outer approximation ---

I = interval([-1; -2; 1], [3; 1; 4]);
E = ellipsoid(I);
resvec(end+1) = contains(E,I);

E = ellipsoid(I, 'outer');
resvec(end+1) = contains(E,I);

% check inner approximation ---

E = ellipsoid(I, 'inner');
resvec(end+1) = contains(I,E);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
