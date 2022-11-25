function res = test_ellipsoid_enlarge
% test_ellipsoid_enlarge - unit test function of enlarge
%
% Syntax:  
%    res = test_ellipsoid_enlarge
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate ellipsoid
E = ellipsoid([1 0;0 2],[0; 1]);
factor = 2;

% compute enlarged ellipsoid
E_enlarged = enlarge(E,factor);

% true solution
E_true = ellipsoid([4 0;0 8],[0; 1]);

% compare results
res = isequal(E_true,E_enlarged);

if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

%------------- END OF CODE --------------