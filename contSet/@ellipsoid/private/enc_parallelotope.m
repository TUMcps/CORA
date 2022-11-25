function Z = enc_parallelotope(E)
% enc_parallelotope - overapproximates an ellipsoid by a parellelotope
%
% Syntax:  
%    Z = enc_parallelotope(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid.generateRandom(0,2);
%    Z = enc_parallelotope(E);
%    plot(E);
%    hold on
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:       Victor Gassmann, Matthias Althoff
% Written:      14-October-2019
% Last update:  27-January-2021 (MA, degenerate case implemented)
%               08-June-2021 (VG, moved degeneracy to main file)
% Last revision:---

%------------- BEGIN CODE --------------
assert(~E.isdegenerate,'Degeneracy should be handled in main file!');

Tinv = sqrtm(E.Q);
n = length(E.Q);
%transform ellipsoid into sphere -> square around sphere -> back transform
Z = zonotope([E.q,Tinv*eye(n)]);



%------------- END OF CODE --------------