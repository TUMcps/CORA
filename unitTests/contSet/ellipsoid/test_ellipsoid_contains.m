function res = test_ellipsoid_contains
% test_ellipsoid_contains - unit test function of contains
%
% Syntax:
%    res = test_ellipsoid_contains
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

res = true;
load cases.mat E_c

% empty set: rewrite using emptySet class
% E_e = ellipsoid.empty(2);
% res = contains(E_c{1}.E1,E_e) && ~contains(E_e,E_c{1}.E1);

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = length(E1.q);
    
    Yd = randPoint(Ed1,2*n);
    if contains(E1,Ed1) && ~contains(E1,Yd)
        res = false;
        break;
    end
    
    if contains(E1,E0) && ~contains(E1,E0)
        res = false;
        break;
    end
    % rewrite using emptySet class
%     E_e = ellipsoid();
%     if contains(E_e,E1) || ~contains(E1,E_e)
%         res = false;
%     end
end


% ellipsoid and zonotope
Q = [13 7; 7 5];
q = [1; 2];
E = ellipsoid(Q,q);

c = [1;1];
G = [1 1 1; 1 -1 0];
Z = zonotope(c,G);

res = res && ~contains(E,Z);

% ------------------------------ END OF CODE ------------------------------
