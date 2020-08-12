function [E] = and(E1,E2,checkint)
% and - overloads & operator to compute the intersects two ellipsoids
%       according to Sec. 2.2.7 in [1] 
%
% Syntax:  
%    [E] = intersect(E1,E2)
%
% Inputs:
%    E1,E2 - Ellipsoid object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1=ellipsoid([1 0; 0 1]);
%    E2=ellipsoid([1 1; 1 1]);
%    E =E1 & E2;
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  15-October-2019
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('checkint','var')
    checkint = true;
end
%We assume that E1,E2 are full-dimensional
if ~isa(E1,'ellipsoid') || ~isa(E2,'ellipsoid')
    error('Wrong input arguments');
end
%check if either ellipsoid is degenerate
if E1.isdegenerate || E2.isdegenerate
    error('E1/E2 degenerate, not supported')
%check if same dimension
elseif length(E1.Q)~=length(E2.Q)
    error('E1/E2 have to have same dimensions');
end

%check if ellipsoids are equal (within tolerance)
if E1==E2
    E = E1;
    return;
end
%For details, see [2]
TOL = min(E1.TOL,E2.TOL);
W1 = inv(E1.Q);
W2 = inv(E2.Q);
q1 = E1.q;
q2 = E2.q;
p = compIntersectionParam(W1,q1,W2,q2);

[~,Q,q] = rootfnc(p,W1,q1,W2,q2);
if any(eig(Q)<0) || (checkint && ~isIntersecting(E1,E2))
    E = ellipsoid;
else
    E = ellipsoid(Q,q);
end

%------------- END OF CODE --------------