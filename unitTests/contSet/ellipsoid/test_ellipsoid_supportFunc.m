function res = test_ellipsoid_supportFunc
% test_ellipsoid_supportFunc - unit test function of supportFunc
%
% Syntax:  
%    res = test_ellipsoid_supportFunc
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    
    res = checkSuppFunc(E1) && checkSuppFunc(Ed1) && checkSuppFunc(E0);
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
end

%-- helper
function res = checkSuppFunc(E)
n = dim(E);
[T,S,~] = svd(E.Q);
s = sqrt(diag(S));
res = true;
for i=1:n
    l = T(:,i);
    [val,x] = supportFunc(E,l);
    ri = abs(val-l'*E.q);
    if ~withinTol(s(i),ri,E.TOL) || ~withinTol(norm(x-E.q),s(i),E.TOL)
        res = false;
        break;
    end
end
end
%------------- END OF CODE --------------
