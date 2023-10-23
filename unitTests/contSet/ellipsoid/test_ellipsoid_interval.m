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
    
    if ~all(contains(interval(E1),Y)) || ~all(contains(interval(Ed1),Yd)) || ~contains(interval(E0),Y0)
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
