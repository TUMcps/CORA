function res = test_ellipsoid_in
% test_ellipsoid_in - unit test function of in
%
% Syntax:  
%    res = test_ellipsoid_in
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
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = length(E1.q);
    
    Yd = randPoint(Ed1,2*n);
    if in(E1,Ed1) && ~in(E1,Yd)
        res = false;
        break;
    end
    
    if in(E1,E0) && ~in(E1,E0)
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
