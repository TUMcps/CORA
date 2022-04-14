function res = test_ellipsoid_zonotope
% test_ellipsoid_zonotope - unit test function of zonotope
%
% Syntax:  
%    res = test_ellipsoid_zonotope
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
    
    Z1 = zonotope(E1,2*n,'i:norm');
    Zd1 = zonotope(Ed1,2*n,'o:norm');
    Z0 = zonotope(E0,2*n,'i:norm:bnd');
    
    if ~in(E1,randPoint(Z1,N,'extreme')) || ...
       ~in(Zd1,randPoint(Ed1,N,'extreme')) || ...
       ~(all(rad(interval(Z0))==0) && all(withinTol(center(Z0),E0.q,E0.TOL)))
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
