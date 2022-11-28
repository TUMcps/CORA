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

load cases.mat E_c  

% empty set
res = supportFunc(ellipsoid(),[1;1],'upper') == -Inf ...
    && supportFunc(ellipsoid(),[1;1],'lower') == Inf;

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    
    res = checkSuppFunc(E1) && checkSuppFunc(Ed1) && checkSuppFunc(E0);
end

end


% Auxiliary functions -----------------------------------------------------
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
