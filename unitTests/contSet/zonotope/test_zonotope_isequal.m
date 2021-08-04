function res = test_zonotope_isequal
% test_zonotope_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_zonotope_isequal
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

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  21-April-2020
%               09-August-2020 (enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------


% 1. Analytical Tests -----------------------------------------------------

% create zonotope
Z1 = zonotope([1, 2, 4;
               5, 6, 0;
              -1, 4, 8]);
Z2 = zonotope([1, 2, 4;
               4, 6, 0;
              -1, 4, 8]);
Z3 = zonotope([1, 2, 0, 4;
               5, 6, 0, 0;
              -1, 4, 0, 8]);

% check result
res_val(1) = isequal(Z1,Z3) && ~isequal(Z1,Z2);

% different order of generators
Z1 = zonotope(ones(2,1),[1 2 5 3 3;
                          2 3 0 4 1]);
Z2 = zonotope(ones(2,1),[2 1 3 5 3;
                          3 2 4 0 1]);

res_val(2) = isequal(Z1,Z2);

% add results
res = all(res_val);

if res
    disp('test_isequal successful');
else
    disp('test_isequal failed');
end

%------------- END OF CODE --------------
