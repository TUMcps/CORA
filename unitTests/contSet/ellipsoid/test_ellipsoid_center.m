function res = test_ellipsoid_center
% test_ellipsoid_center - unit test function of center
%
% Syntax:  
%    res = test_ellipsoid_center
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
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate ellipsoids
E1 = ellipsoid([1 0;0 2],[0; 1]);
E2 = ellipsoid([1 0;0 2]); % center at origin

% compare results
res = isequal(center(E1),[0;1]) && isequal(center(E2),[0;0]);

if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

%------------- END OF CODE --------------