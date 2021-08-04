function res = test_ellipsoid_plus
% test_ellipsoid_plus - unit test function of plus
%
% Syntax:  
%    res = test_ellipsoid_plus
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
    N = 2*dim(E1);
    
    Y = randPoint(E1,N);
    Yd = randPoint(Ed1,N);
    Y0 = randPoint(E0,N);
    
    if ~in(E1+Ed1,Y+Yd) || ~in(Ed1+E0,Yd+Y0)
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
