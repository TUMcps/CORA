function res = test_ellipsoid_radius
% test_ellipsoid_radius - unit test function of radius
%
% Syntax:  
%    res = test_ellipsoid_radius
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = dim(E1);
    
    if ~all(withinTol([radius(E0);radius(E0,randi([1,n]))],0,E0.TOL))
        res = false;
        break;
    end
    
    if ~withinTol(norm(ellipsoid(Ed1.Q)),radius(Ed1),Ed1.TOL)
        res = false;
        break;
    end
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
