function res = test_ellipsoid_minus
% test_ellipsoid_minus - unit test function of minus
%
% Syntax:  
%    res = test_ellipsoid_minus
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
    Ed1 = E_c{i}.Ed1;
    E0 = E_c{i}.E0;
    
    % since E0.Q = 0
    if isempty(minus(E1,E0))
        res = false;
        break;
    end
    
    Ep = E0+Ed1;
    if ~in((Ep-Ed1)+Ed1,Ep)
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
