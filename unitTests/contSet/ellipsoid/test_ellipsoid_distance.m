function res = test_ellipsoid_distance
% test_ellipsoid_distance - unit test function of distance
%
% Syntax:
%    res = test_ellipsoid_distance
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
% Last update:   07-July-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    E0 = E_c{i}.E0; % all zero
    
    % check all-zero ellipsoid
    if ~withinTol(distance(E1,E0),distance(E1,E0.q),E1.TOL)
        res = false;
        break;
    end
    
    n = length(E1.q);
    % check conHyperplane: construct hyperplane through E1.q
    l1 = randn(n,1);
    H = conHyperplane(l1,l1'*E1.q);
    % check if distance is <0
    if distance(E1,H)>=0
        res = false;
        break;
    end
    
    % check polytope: construct second hyperplane (also contains center
    % of ellipsoid)
    l2 = randn(n,1);
    P = polytope([l1';l2'],[l1'*E1.q;l2'*E1.q]);
    if distance(E1,P)>1e-6
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
