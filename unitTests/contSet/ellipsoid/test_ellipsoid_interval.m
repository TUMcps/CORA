function res = test_ellipsoid_interval
% test_ellipsoid_interval - unit test function of interval
%
% Syntax:  
%    res = test_ellipsoid_interval
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
    n = length(E1.q);
    
    Y = randPoint(E1,2*n);
    Yd = randPoint(Ed1,2*n);
    Y0 = E0.q;
    
    if ~in(interval(E1),Y) || ~in(interval(Ed1),Yd) || ~in(interval(E0),Y0)
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
