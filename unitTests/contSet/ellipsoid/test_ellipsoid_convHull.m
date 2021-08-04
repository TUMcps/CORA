function res = test_ellipsoid_convHull
% test_ellipsoid_convHull - unit test function of convHull
%
% Syntax:  
%    res = test_ellipsoid_convHull
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

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
[~,res] = evalc('test_ellipsoid_or');


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
