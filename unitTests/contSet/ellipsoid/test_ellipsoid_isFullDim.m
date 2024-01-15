function res = test_ellipsoid_isFullDim
% test_ellipsoid_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_ellipsoid_isFullDim
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

load cases.mat E_c

% empty set
res = ~isFullDim(ellipsoid.empty(2));

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = length(E1.q);
    
    if ~isFullDim(E1) || isFullDim(Ed1) || isFullDim(E0)
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
