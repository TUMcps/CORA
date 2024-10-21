function [res,val] = priv_containsPoint(E,S)
% priv_containsPoint - checks whether an ellipsoid contains a point cloud;
%    additionally, the relative distance to the center of the ellipsoid
%    is returned (<= 1 if contained, =1 on boundary, otherwise Inf)
%
% Syntax:
%    [res,val] = priv_containsPoint(E,S)
%
% Inputs:
%    E1 - ellipsoid object
%    S - point cloud
%
% Outputs:
%    res - true/false for each point
%    val - distance
%
% References:
%    [1] Boyd et al. Convex Optimization (B.2, ex. B.1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/contains_

% Authors:       Victor Gassmann
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of points
numPoints = size(S,2);
% logical mask for those that still need to be checked
mask = true(1,numPoints);

% init results
res = false(1,numPoints);
val = Inf(1,numPoints);

% degenerate case
if ~isFullDim(E)

    % subspace dimension
    n_subspace = rank(E);

    % all-zero Q matrix -> ellipsoid is just a point, so we just check for
    % equality between the points and the center (distance is 0 or Inf)
    if n_subspace == 0
        res = all(withinTol(E.q,S,E.TOL),1);
        val(res) = 0;
        return;
    end

    % project ellipsoid and points onto subspace via SVD
    [T,~,~] = svd(E.Q);
    E = T'*E;
    S = T'*S;

    % save remaining flat dimensions
    q_remainder = E.q(n_subspace+1:end);
    S_remainder = S(n_subspace+1:end,:);

    % check whether flat dimensions have same value (necessary condition
    % for any point to be contained in the ellipsoid)
    ind_rem_eq = all(withinTol(q_remainder,S_remainder,E.TOL),1);
    mask = mask & ind_rem_eq;

    % if only center remains, we can exit immediately
    if rank(E)==0
        res(ind_rem_eq) = true;
        val(ind_rem_eq) = 1;
        return;
    end

    % project onto first few dimensions so that E is no longer degenerate
    rankE = rank(E);
    E = project(E,1:rankE);
    S = S(1:rankE,:);
end

% now, E is full-dimensional
tol = 1 + E.TOL;
for i=1:numPoints
    if mask(i)
        % simply check using ellipsoid equation
        val_i = (S(:,i)-E.q)' * (E.Q \ (S(:,i)-E.q));
        res(i) = val_i < tol | withinTol(val_i,tol);
        val(i) = res(i)*val_i;
    end
end

% ------------------------------ END OF CODE ------------------------------
