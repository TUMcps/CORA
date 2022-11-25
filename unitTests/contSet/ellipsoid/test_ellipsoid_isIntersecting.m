function res = test_ellipsoid_isIntersecting
% test_ellipsoid_isIntersecting - unit test function of isIntersecting
%
% Syntax:  
%    res = test_ellipsoid_isIntersecting
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
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.E1;
    E0 = E_c{i}.E0;
    
    E0 = ellipsoid(E0.Q,E1.q);
    
    if ~isIntersecting(E1,E0)
        res = false;
        break;
    end
    
    r1 = rad(interval(E1));
    rd1 = rad(interval(Ed1));
    Ed1 = ellipsoid(Ed1.Q,Ed1.q+r1+rd1);
    if isIntersecting(Ed1,E1)
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
