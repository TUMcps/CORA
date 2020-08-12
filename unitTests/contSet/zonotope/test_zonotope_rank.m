function res = test_zonotope_rank
% test_zonotope_rank - unit test function of rank
%
% Syntax:  
%    res = test_zonotope_rank
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
% Last update:  15-Sep-2019
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z1 = zonotope([1, 2, 0, 4; 5, 6, 0, 0; -1, 4, 0, 8]);

% obtain zonotope without zeros
d = rank(Z1);

% true result
true_val = 2;

% check result
res = (d==true_val);


if res
    disp('test_rank successful');
else
    disp('test_rank failed');
end

%------------- END OF CODE --------------
