function res = test_ellipsoid_randPoint
% test_ellipsoid_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_ellipsoid_randPoint
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
    N = 5*n;
    
    if ~withinRange(E1,N) || ~withinRange(Ed1,N) || ~withinRange(E0,N)
        res = false;
        break;
    end
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end

end

%-- helper
function res = withinRange(E,N)
    % all extreme points need to be between min(radius) and max(radius)
    n = dim(E);
    Y = randPoint(E,N,'extreme');
    nY = sqrt(sum((Y-E.q).^2,1));
    rE = radius(E,n);
    IntR = interval(min(rE),max(rE));
    res = in(IntR,nY);
end
%------------- END OF CODE --------------
