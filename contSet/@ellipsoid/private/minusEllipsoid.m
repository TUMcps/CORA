function E = minusEllipsoid(E1,E2,L,mode)
% minusEllipsoid - Computes the inner or outer approximation of the
% difference between two ellipsoids
%
% Syntax:  
%    E = minusEllipsoid(E,E2,L,mode)
%
% Inputs:
%    E1,E2  - ellipsoid objects
%    L      - directions 
%    mode   - mode ('i':inner approx; 'o': outer approx)
%
% Outputs:
%    E - ellipsoid after Minkowski difference
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      15-March-2021
% Last update:  27-July-2021 (VG: added zero cases)
% Last revision:---

%------------- BEGIN CODE --------------
if rank(E1)==0 && rank(E2)==0
    E = ellipsoid(zeros(dim(E1)),E1.q-E2.q);
    return;
elseif rank(E1)==0 && rank(E2)>0
    E = ellipsoid;
    return;
elseif rank(E1)>0 && rank(E2)==0
    E = ellipsoid(E1.Q,E1.q-E2.q);
    return;
end

% check if difference should be empty
if ~in(ellipsoid(E1.Q),ellipsoid(E2.Q))
    E = ellipsoid;
    return;
end
n = E1.dim;
TOL = min(E1.TOL,E2.TOL);
if all(all(withinTol(E1.Q,E2.Q,TOL)))
    E = ellipsoid(zeros(n),E1.q-E2.q);
    return;
end

% check bad directions
isbDir = @(L) isBadDir(L,E1,E2);
if isempty(L)
    N = 2*n;
    if n==1
        L = [-1,1];
        L(:,isbDir(L)) = [];
    else
        L = eq_point_set(n-1,N);
        counter = 1;
        while size(L,2)<N
            L = eq_point_set(n-1,N+counter);
            L(:,isbDir(L)) = [];
            counter = counter + 1;
        end
        L = L(:,1:N);
    end
    
else
    L(:,isbDir(L)) = [];
    if isempty(L)
        error('Provided directions are all not suitable for the Minkowski difference');
    end
end
E_cell = lminus(E1,E2,L,'o');
if strcmp(mode','o')
    % compute intersection (outer approx)
    E = E_cell{1};
    E = and(E,E_cell(2:end),'o');
else
    % bulk inner union
    E = or(E_cell{1},E_cell(2:end),'i');
end
%------------- END OF CODE --------------