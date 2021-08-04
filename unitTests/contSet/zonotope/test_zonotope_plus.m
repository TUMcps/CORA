function res = test_zonotope_plus
% test_zonotope_plus - unit test function of plus
%
% Syntax:  
%    res = test_zonotope_plus
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical Test ------------------------------------------------------

% create zonotopes
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);
Z2 = zonotope([1 10; -1 -10]);

% obtain results
Z3 = Z1+Z2;

% obtain zonotope matrix
Zmat = Z3.Z;

% true result
true_mat = [-3, -3, -2, -1, 10; ...
            0, 2, 3, 4, -10];

% check result
res_val = all(all(Zmat == true_mat));

% add results
res = res_val;

if res
    disp('test_zonotope_plus successful');
else
    disp('test_zonotope_plus failed');
end

%------------- END OF CODE --------------
