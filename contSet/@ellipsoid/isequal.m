function res = isequal(E1,E2)
% isequal - checks if ellipsoids E1, E2 are equal within tolerance
%
% Syntax:  
%    B = isequal(E,Y) gives a logical indicating whether E1,E2 are equal
%    (within tolerance)
%
% Inputs:
%    E1 - ellipsoids object
%    E2 - ellipsoids object
%
% Outputs:
%    res - boolean value indicating whether points E1,E2 are equal
%
% Example: 
%    E1 = ellipsoid([1,0;0,1/2],[1;1]);
%    E2 = ellipsoid([1+1e-15,0;0,1/2],[1;1]);
%    B = equal(E1,E2);
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
TOL = min(E1.TOL,E2.TOL);
res_q = all(abs(E1.q-E2.q)<=TOL);
res_Q = all(all(abs(E1.Q-E2.Q)<=TOL));
res = res_q && res_Q;

%------------- END OF CODE --------------