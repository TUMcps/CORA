function res = test_zonotope_deleteZeros
% test_deleteZeros - unit test function of deleteZeros
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([1,2,0,4; 5 6 0 0]);

% obtain zonotope without zeros
Z2 = deleteZeros(Z1);

% obtain zonotope matrix
Zmat = get(Z2,'Z');

% true result
true_mat = [1, 2, 4; 5, 6, 0];

% check result
res = all(all(Zmat == true_mat));


if res
    disp('test_deleteZeros successful');
else
    disp('test_deleteZeros failed');
end

%------------- END OF CODE --------------
