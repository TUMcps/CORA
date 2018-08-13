function res = test_zonotope_center
% test_zonotope_center - unit test function of center
%
% Syntax:  
%    res = test_center
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
Z1 = zonotope([1,2,3,4; 5 6 7 8]);

% obtain center
c = center(Z1);

% true result
true_vec = [1; 5];

% check result
res = all(c == true_vec);

if res
    disp('test_center successful');
else
    disp('test_center failed');
end

%------------- END OF CODE --------------
