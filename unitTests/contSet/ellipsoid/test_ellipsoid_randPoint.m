function res = test_ellipsoid_randPoint
% test_ellipsoid_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_ellipsoid_randPoint
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
n = 2;
E = ellipsoid.empty(n);
p = randPoint(E);
res = isempty(p) && isnumeric(p) && size(p,1) == n;

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = dim(E1);
    N = 5*n;
    
    if ~aux_withinRange(E1,N) || ~aux_withinRange(Ed1,N) || ~aux_withinRange(E0,N)
        res = false;
        break;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_withinRange(E,N)
    % all extreme points need to be between min(radius) and max(radius)
    n = dim(E);
    Y = randPoint(E,N,'extreme');
    nY = sqrt(sum((Y-E.q).^2,1));
    rE = radius(E,n);
    IntR = interval(min(rE),max(rE));
    res = all(contains(IntR,nY));
end

% ------------------------------ END OF CODE ------------------------------
