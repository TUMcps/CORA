function res = test_ellipsoid_or
% test_ellipsoid_or - unit test function of or
%
% Syntax:
%    res = test_ellipsoid_or
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

% empty set
E_empty = ellipsoid.empty(2);
if ~isequal(or(E_c{1}.E1,E_empty),E_c{1}.E1)
    res = false;
end

for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    E2 = E_c{i}.E2; % non-deg
    E0 = E_c{i}.E0; % all zero
    
    % test non-deg
    Eres_nd = E1 | E2;
    Y_nd = [randPoint(E1,2*i),randPoint(E2,2*i)];
    if ~contains(Eres_nd,Y_nd) 
        res = false;
        break;
    end
    
    % test zero rank ellipsoid
    Eres_0 = E1 | E0;
    Y_0 = [randPoint(E1,2*i),E0.q];
    if ~contains(Eres_0,Y_0)
        res = false;
        break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
