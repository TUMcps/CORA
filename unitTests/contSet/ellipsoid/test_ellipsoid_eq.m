function res = test_ellipsoid_eq
% test_ellipsoid_eq - unit test function of eq
%
% Syntax:
%    res = test_ellipsoid_eq
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

% Authors:       Victor Gassmann
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[~,res] = evalc('test_ellipsoid_isequal');

% ------------------------------ END OF CODE ------------------------------
