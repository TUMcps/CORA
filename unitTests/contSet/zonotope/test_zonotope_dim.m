function res = test_zonotope_dim
% test_dim - unit test function of dim
%
% Syntax:  
%    res = test_zonotope_dim
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
Z1 = zonotope([1, 2, 0, 4; 5, 6, 0, 0; -1, 4, 0, 8]);

% obtain zonotope without zeros
d = dim(Z1);

% true result
true_val = 2;

% check result
res = (d==true_val);


if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------
