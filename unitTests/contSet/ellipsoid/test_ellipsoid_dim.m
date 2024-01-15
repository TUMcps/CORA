function res = test_ellipsoid_dim
% test_ellipsoid_dim - unit test function of dim
%
% Syntax:
%    res = test_ellipsoid_dim
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
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty set
n = 2;
E = ellipsoid.empty(n);
res = dim(E) == n;

load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    
    if ~all([dim(E1),dim(Ed1),dim(E0)]==length(E1.q))
        res = false;
        break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
