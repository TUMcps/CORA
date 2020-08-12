function Z = insc_parallelotope(E)
% insc_parallelotope - underapproximates an ellipsoid by a parellelotope
%
% Syntax:  
%    Z = insc_parallelotope(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid.generateRandom(0,2);
%    Z = insc_parallelotope(E);
%    plot(E);
%    hold on
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:       Victor Gassmann
% Written:      14-October-2019
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
if E.isdegenerate
    error('not implemented for degenerate ellipsoids');
end
T = inv(sqrtm(E.Q));
n = length(E.Q);
%transform ellipsoid into sphere -> square into sphere -> back transform
Z = zonotope([E.q,inv(T)*1/sqrt(n)*eye(n)]);
%------------- END OF CODE --------------