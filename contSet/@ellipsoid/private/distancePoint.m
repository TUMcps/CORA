function D = distancePoint(E,Y)
% distancePoint - computes the distance from an ellipsoid to an (array of)
%    points
%
% Syntax:
%    D = distancePoint(E,Y)
%
% Inputs:
%    E - ellipsoid object
%    Y - point or matrix of points
%
% Outputs:
%    D - distance(s) between E and Y
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   19-May-2022
%                02-June-2022 (VG, complete rework)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% say Q is given by blkdiag(Q,zeros). We look at is as blkdiag(Q,TOL).
% => (x-q)'*Q\(x-q) + (y_rem-x_rem)'*1/TOL*(y_rem-x_rem) <= 1
N = size(Y,2);
D_d = zeros(1,N);
if ~isFullDim(E)
    nt = rank(E);
    if nt==0
        % if only center remains, compute distance using ellipsoid equation
        D = sum((1/sqrt(E.TOL)*(Y-E.q)).^2,1) - 1;
        return;
    end
    [T,~,~] = svd(E.Q);    
    % transform E and Y (so that degeneracy of E is axis-aligned)
    E = T'*E;
    Y = T'*Y;
    Y_rem = Y(nt+1:end,:);
    x_rem = E.q(nt+1:end);
    E = project(E,1:nt);
    Y = Y(1:nt,:);
    % precompute distance from x_rem
    D_d = sum((1/sqrt(E.TOL)*(Y_rem-x_rem)).^2,1);
end

D_nd = sum((sqrtm(E.Q)\(Y-E.q)).^2,1);
% concatenate and subtract 1 (so that 0 means touching (instead of 1))
D = D_nd + D_d - 1;

% ------------------------------ END OF CODE ------------------------------
