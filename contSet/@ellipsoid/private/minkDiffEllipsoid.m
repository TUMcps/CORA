function E = minkDiffEllipsoid(E1,E2,L,mode)
% minkDiffEllipsoid - Computes an inner- or outer-approximation of the
%    Minkowski difference between two ellipsoids
%
% Syntax:
%    E = minkDiffEllipsoid(E,E2,L,mode)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object
%    L - directions 
%    mode - type of approximation: 'inner', 'outer'
%
% Outputs:
%    E - ellipsoid after Minkowski difference
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference on
%       Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/minkDiff

% Authors:       Victor Gassmann
% Written:       15-March-2021
% Last update:   27-July-2021 (VG, added zero cases)
%                09-November-2022 (MW, rename 'minkDiffEllipsoid')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simply center difference if both are simply their respective centers
if rank(E1)==0 && rank(E2)==0
    E = ellipsoid(zeros(dim(E1)),E1.q-E2.q);
    return;
% since we check if E2\subseteq E1, this cannot hold if E2 has higher rank
elseif rank(E1)==0 && rank(E2)>0
    E = ellipsoid;
    return;
% simply subtract E2 (which is a point)
elseif rank(E1)>0 && rank(E2)==0
    E = ellipsoid(E1.Q,E1.q-E2.q);
    return;
end

% check if difference should be empty (only non-empty if
% ellipsoid(E2.Q)\subseteq ellipsoid(E1.Q)
if ~isBigger(E1,E2)
    E = ellipsoid.empty(dim(E1));
    return;
end
n = dim(E1);
TOL = min(E1.TOL,E2.TOL);
% check if shape matrices match (within tolerance)
if all(all(withinTol(E1.Q,E2.Q,TOL)))
    E = ellipsoid(zeros(n),E1.q-E2.q);
    return;
end

% check bad directions and remove them
isbDir = @(L) isBadDir(E1,E2,L);
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
        throw(CORAerror('CORA:emptySet'));
%         error('Provided directions are all not suitable for the Minkowski difference');
    end
end
E_cell = lminkDiff(E1,E2,L,mode);
if length(E_cell)==1
    E = E_cell{1};
elseif strcmp(mode','outer')
    % compute intersection (outer approx)
    E = E_cell{1};
    E = and_(E,[E_cell{2:end}],'outer');
else
    % bulk inner union
    E = or(E_cell{1},[E_cell{2:end}],'inner');
end

% ------------------------------ END OF CODE ------------------------------
