function res = test_ellipsoid_volume
% test_ellipsoid_volume - unit test function of volume
%
% Syntax:
%    res = test_ellipsoid_volume
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% res = true;
res = (volume(ellipsoid.empty(2)) == 0);
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = dim(E1);
    
    E1_vol = pi^(n/2)/gamma(n/2+1)*sqrt(det(E1.Q));
    
    if ~withinTol(volume(E0),0,E0.TOL) || ~withinTol(volume(Ed1),0,Ed1.TOL) ||...
            ~withinTol(volume(E1),E1_vol,E1.TOL)
        res = false;
        break;
    end    
end

% ------------------------------ END OF CODE ------------------------------
