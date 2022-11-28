function Z = insc_parallelotope(E)
% insc_parallelotope - inner-approximates an ellipsoid by a parellelotope
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
%    E = ellipsoid.generateRandom('Dimension',2);
%    Z = insc_parallelotope(E);
%
%    figure; hold on;
%    plot(E);
%    plot(Z,[1,2],'r');
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

if ~isFullDim(E)
    throw(CORAerror('CORA:degenerateSet','Should be handled in main file'));
end

T = inv(sqrtm(E.Q));
n = length(E.Q);
%transform ellipsoid into sphere -> square into sphere -> back transform
Z = zonotope([E.q,inv(T)*1/sqrt(n)*eye(n)]);

%------------- END OF CODE --------------