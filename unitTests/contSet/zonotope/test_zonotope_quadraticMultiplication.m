function res = test_zonotope_quadraticMultiplication
% test_zonotope_quadraticMultiplication - unit test function of quadraticMultiplication
%
% Syntax:  
%    res = test_zonotope_quadraticMultiplication
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
Z1 = zonotope([-4, -3, -2; 1, 2, 3]);

% create matrices
Q{1} = [1 -2; -3 4];
Q{2} = [0.5 0; 2 -1];

% obtain result
Z2 = quadraticMultiplication(Z1,Q);

% obtain zonotope matrix
Zmat = get(Z2,'Z');

% true result
true_mat = [102.5, 27.5, 35, 95, 110, 125; ...
            -16.25, -5.75, -9.5, -14, -26, -32];

% check result
res = all(all(Zmat == true_mat));

if res
    disp('test_zonotope_quadraticMultiplication successful');
else
    disp('test_zonotope_quadraticMultiplication failed');
end

%------------- END OF CODE --------------
