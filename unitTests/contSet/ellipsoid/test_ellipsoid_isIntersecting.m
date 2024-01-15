function res = test_ellipsoid_isIntersecting
% test_ellipsoid_isIntersecting - unit test function of isIntersecting
%
% Syntax:
%    res = test_ellipsoid_isIntersecting
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
res = true;

% empty set: rewrite using emptySet class
% res = ~ isIntersecting(E_c{1}.E1,ellipsoid.empty(2));

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.E1;
    E0 = E_c{i}.E0;
    
    E0 = ellipsoid(E0.Q,E1.q);
    
    if ~isIntersecting(E1,E0)
        res = false;
        break;
    end
    
    r1 = rad(interval(E1));
    rd1 = rad(interval(Ed1));
    Ed1 = ellipsoid(Ed1.Q,Ed1.q+r1+rd1);
    if isIntersecting(Ed1,E1)
        res = false;
        break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
