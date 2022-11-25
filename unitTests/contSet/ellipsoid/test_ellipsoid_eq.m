function res = test_ellipsoid_eq
% test_ellipsoid_eq - unit test function of eq
%
% Syntax:  
%    res = test_ellipsoid_eq
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
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

[~,res] = evalc('test_ellipsoid_isequal');


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
