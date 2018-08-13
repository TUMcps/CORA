function res = test_zonotope_project
% test_zonotope_project - unit test function of project
%
% Syntax:  
%    res = test_zonotope_project
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
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4; 5, 5, 5, 5]);

% obtain result
Z2 = project(Z1,[1 3]);

% obtain zonotope matrix
Zmat = get(Z2,'Z');

% true result
true_mat = [-4, -3, -2, -1; ...
            5, 5, 5, 5];

% check result
res = all(all(Zmat == true_mat));

if res
    disp('test_zonotope_project successful');
else
    disp('test_zonotope_project failed');
end

%------------- END OF CODE --------------
