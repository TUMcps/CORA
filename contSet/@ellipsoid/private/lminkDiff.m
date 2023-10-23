function E_cell = lminkDiff(E1,E2,L,mode)
% lminkDiff - Approximate the Minkowski difference of an ellipsoid and
%    another ellipsoid or a vector
%
% Syntax:
%    E_cell = lminkDiff(E1,E2,L,mode)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object or numerical vector
%    L - unit vectors in different directions
%    mode - type of approximation: 'inner', 'outer'
%
% Outputs:
%    E_cell - cell array of ellipsoid after Minkowski difference (length ==
%               size(L,2))
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
% See also: -

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   13-June-2022 (VG, Bugfix)
%                09-November-2022 (MW, rename 'lminkDiff')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if we have a bad direction in L
if any(isBadDir(E1,E2,L))
    throw(CORAerror('CORA:wrongValue','third',...
        "'L' should not contain bad direction (see 'isBadDir.m')"));
end


% make sure L is normalized
L = 1./sqrt(sum(L.^2,1)).*L;

E_cell = cell(1,size(L,2));
% continue with "good" directions
for i=1:size(L,2)
    q = E1.q-E2.q;
    l = L(:,i);
    if strcmp(mode,'outer')
        Q1s = sqrtm(E1.Q);
        Q2s = sqrtm(E2.Q);
        % there is a sign error in the pdf manual of the ET
        Q_ = Q1s-vecalign(Q1s*l,Q2s*l)*Q2s;
        Q = Q_'*Q_;
    else
        lql1 = sqrt(l'*E1.Q*l);
        lql2 = sqrt(l'*E2.Q*l);
        p = lql1/lql2;
        Q = (1-1/p)*E1.Q + (1-p)*E2.Q;
    end
    E_cell{i} = ellipsoid(Q,q);
end

% ------------------------------ END OF CODE ------------------------------
