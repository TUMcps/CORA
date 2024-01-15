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
E = ellipsoid.empty(2);
res = supportFunc(E,[1;1],'upper') == -Inf ...
    && supportFunc(E,[1;1],'lower') == Inf;

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    
    res = aux_checkSuppFunc(E1) && aux_checkSuppFunc(Ed1) && aux_checkSuppFunc(E0);
end

% check type = 'range'
E = ellipsoid([5 7;7 13],[1;2]);
dir = [1;1];

[val_upper,x_upper] = supportFunc(E,dir);
[val_lower,x_lower] = supportFunc(E,dir,'lower');
[val_int,x_both] = supportFunc(E,dir,'range');
if ~isequal(interval(val_lower,val_upper),val_int) ...
        || ~compareMatrices([x_upper, x_lower],x_both)
    res = false;
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkSuppFunc(E)
    n = dim(E);
    [T,S,~] = svd(E.Q);
    s = sqrt(diag(S));
    res = true;

    % loop over all directions
    for i=1:n
        % direction
        l = T(:,i);
        % evaluate support function and compute support vector
        [val,x] = supportFunc(E,l);
        ri = abs(val-l'*E.q);

        % check results
        if ~withinTol(s(i),ri,E.TOL) || ~withinTol(norm(x-E.q),s(i),E.TOL)
            res = false;
            break;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
