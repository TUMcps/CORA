function res = test_zonotope_enlarge
% test_enlarge - unit test function of enlarge
%
% Syntax:  
%    res = test_zonotope_enlarge
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

% obtain enlarged zonotope
Z2 = enlarge(Z1,[2,1.5]);

% obtain zonotope matrix
Zmat = get(Z2,'Z');

% true result
true_mat = [-4, -6, -4, -2; ...
            1, 3, 4.5, 6];

% check result
res = all(all(Zmat == true_mat));

if res
    disp('test_enlarge successful');
else
    disp('test_enlarge failed');
end

%------------- END OF CODE --------------
