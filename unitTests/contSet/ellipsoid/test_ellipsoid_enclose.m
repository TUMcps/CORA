function res = test_ellipsoid_enclose
% test_ellipsoid_enclose - unit test function of enclose
%
% Syntax:  
%    res = test_ellipsoid_enclose
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
    E2 = E_c{i}.E2; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = length(E1.q);
    
    % test non-deg ellipsoid
    E = enclose(E1,E2);
    Y = [randPoint(E1,2*n),randPoint(E2,2*n)];
    if ~in(E,Y)
        res = false;
        break;
    end
    
    % test degenerate ellipsoids
    Ed = enclose(Ed1,E0);
    Yd = [randPoint(Ed1,2*n),E0.q];
    if ~in(Ed,Yd)
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
