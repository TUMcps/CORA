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

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);
Ed1 = ellipsoid([ 4.2533342807136076 0.6346400221575308 ; 0.6346400221575309 0.0946946398147988 ], [ -2.4653656883489115 ; 0.2717868749873985 ], 0.000001);
E0 = ellipsoid([ 0.0000000000000000 0.0000000000000000 ; 0.0000000000000000 0.0000000000000000 ], [ 1.0986933635979599 ; -1.9884387759871638 ], 0.000001);

% empty set
E = ellipsoid.empty(2);
assert(supportFunc(E,[1;1],'upper') == -Inf ...
    && supportFunc(E,[1;1],'lower') == Inf);
    
assert(aux_checkSuppFunc(E1))
assert(aux_checkSuppFunc(Ed1))
assert(aux_checkSuppFunc(E0));

% check type = 'range'
E = ellipsoid([5 7;7 13],[1;2]);
dir = [1;1];

[val_upper,x_upper] = supportFunc(E,dir);
[val_lower,x_lower] = supportFunc(E,dir,'lower');
[val_int,x_both] = supportFunc(E,dir,'range');
assert(isequal(interval(val_lower,val_upper),val_int))
assert(compareMatrices([x_upper, x_lower],x_both))

res = true;

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
        assert(withinTol(s(i),ri,E.TOL))
        assert(withinTol(norm(x-E.q),s(i),E.TOL))
    end
end

% ------------------------------ END OF CODE ------------------------------
