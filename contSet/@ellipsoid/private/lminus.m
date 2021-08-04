function E_cell = lminus(E1,E2,L,mode)
% minus - Approximate the Minkowski difference of two
% ellipsoids 
%
% Syntax:  
%    E = minus(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object or numerical vector
%    L  - unit vectors in different directions
%
% Outputs:
%    E - ellipsoid after Minkowski difference
%
% Example: 
%
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      09-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('mode','var')
    mode = 'o';
end
if ~strcmp(mode,'o') && ~strcmp(mode,'i')
    error('Mode has to be either "o" (outer) or "i" (inner)!');
end

% check if we have a bad direction in L
if any(isBadDir(L,E1,E2))
    error('"L" contains at least one bad direction (see "isBadDir.m")');
end


% make sure L is normalized
L = 1./sqrt(sum(L.^2,1)).*L;

E_cell = cell(1,size(L,2));
% continue with "good" directions
for i=1:size(L,2)
    q = E1.q-E2.q;
    l = L(:,i);
    if strcmp(mode,'o')
        Q1s = sqrtm(E1.Q);
        Q2s = sqrtm(E2.Q);
        Q_ = Q1s+vecalign(Q1s*l,Q2s*l)*Q2s;
        Q = Q_'*Q_;
    else
        lql1 = sqrt(l'*E1.Q*l);
        lql2 = sqrt(l'*E2.Q*l);
        Q = (1-lql1/lql2)*E1.Q + (1-lql2/lql1)*E2.Q;
    end
    E_cell{i} = ellipsoid(Q,q);
end
%------------- END OF CODE --------------