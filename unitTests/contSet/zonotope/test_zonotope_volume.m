function res = test_zonotope_volume
% test_zonotope_volume - unit test function of volume
%
% Syntax:  
%    res = test_zonotope_volume
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

% obtain result
vol = volume(Z1);

% true result 1
true_vol = 80;
        

% check results
res = (vol == true_vol);


if res
    disp('test_zonotope_volume successful');
else
    disp('test_zonotope_volume failed');
end

%------------- END OF CODE --------------
