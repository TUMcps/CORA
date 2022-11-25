function res = test_zonotope_center
% test_zonotope_center - unit test function of center
%
% Syntax:  
%    res = test_zonotope_center
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
% Last update:  09-August-2020 (MW, enhance randomness of test)
% Last revision:---

%------------- BEGIN CODE --------------


% 1. Analytical Test ------------------------------------------------------

% create zonotope
Z = zonotope([1,2,3,4; 5 6 7 8]);

% obtain center
c = center(Z);

% true result
true_vec = [1; 5];

% check result
res_val = all(c == true_vec);

% add results
res = res_val;

if res
    disp('test_center successful');
else
    disp('test_center failed');
end

%------------- END OF CODE --------------
