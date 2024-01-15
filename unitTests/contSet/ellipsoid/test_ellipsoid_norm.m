function res = test_ellipsoid_norm
% test_ellipsoid_norm - unit test function of norm
%
% Syntax:
%    res = test_ellipsoid_norm
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
res = norm(ellipsoid.empty(2)) == -Inf;

% loop over cases
for i=1:length(E_c)
    E1 = ellipsoid(E_c{i}.E1.Q);
    Ed1 = ellipsoid(E_c{i}.Ed1.Q);
    E0 = ellipsoid(E_c{i}.E0.Q);
    TOL = E1.TOL;
    n = dim(E1);
    N = 10*n;
    norm1 = max(sqrt(sum(randPoint(E1,N,'extreme').^2,1)));
    normd = max(sqrt(sum(randPoint(Ed1,N,'extreme').^2,1)));
    norm0 = max(sqrt(sum(randPoint(E0,N,'extreme').^2,1)));
    
    if norm(E1)+TOL<norm1 || norm(Ed1)+TOL<normd || norm(E0)+TOL<norm0
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
