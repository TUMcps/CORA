function res = test_zonotope_mtimes
% test_zonotope_mtimes - unit test function of mtimes
%
% Syntax:  
%    res = test_zonotope_mtimes
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% create parallelotopes
M = [-1 2; 3 -4];

% obtain results
Z2 = M*Z1;

% obtain zonotope matrix
Zmat = get(Z2,'Z');

% true result
true_mat = [6, 7, 8, 9; ...
            -16, -17, -18, -19];

% check result
res = all(all(Zmat == true_mat));

if res
    disp('test_zonotope_mtimes successful');
else
    disp('test_zonotope_mtimes failed');
end

%------------- END OF CODE --------------
