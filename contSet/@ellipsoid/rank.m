function r = rank(E)
% rank - returns the rank of E
%
% Syntax:  
%    d = rank(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    r - rank of E
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    r = rank(E) 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% Find minimum svd threshold using reciprocal
% condition number
d = svd(E.Q);
mev_th = max(d)*E.TOL;
r = sum(d>0 & d>=mev_th);
%------------- END OF CODE --------------