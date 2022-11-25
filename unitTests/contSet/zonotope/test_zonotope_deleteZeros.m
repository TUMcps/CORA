function res = test_zonotope_deleteZeros
% test_zonotope_deleteZeros - unit test function of deleteZeros
%    this encompasses checking the function nonzeroFilter
%
% Syntax:  
%    res = test_zonotope_deleteZeros
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
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Analytical Test ------------------------------------------------------

% create zonotope
Z1 = zonotope([1,2,0,4; 5 6 0 0]);

% obtain zonotope without zeros
Z2 = deleteZeros(Z1);

% obtain zonotope matrix
Zmat = Z2.Z;

% true result
true_mat = [1, 2, 4; 5, 6, 0];

% check result
res_val = all(all(Zmat == true_mat));




% add results
res = res_val;

if res
    disp('test_deleteZeros successful');
else
    disp('test_deleteZeros failed');
end

%------------- END OF CODE --------------
