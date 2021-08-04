function res = test_zonotope_cartProd
% test_zonotope_cartProd - unit test function of cartesian product
%
% Syntax:  
%    res = test_zonotope_cartProd
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

% create zonotopes
Z1 = zonotope([1,2,3,4; 5 6 7 8]);
Z2 = zonotope([9 10 11]);

% compute Cartesian product
Z3 = cartProd(Z1,Z2);

% obtain zonotope matrix
Zmat = Z3.Z;

% true result
true_mat = [1, 2, 3, 4, 0, 0; ...
            5, 6, 7, 8, 0, 0; ...
            9, 0, 0, 0, 10,11];

% check result
res = all(all(Zmat == true_mat));

if res
    disp('test_cartProd successful');
else
    disp('test_cartProd failed');
end

%------------- END OF CODE --------------
