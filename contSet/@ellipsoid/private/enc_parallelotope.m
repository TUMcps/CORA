function Z = enc_parallelotope(E)
% enc_parallelotope - over-approximates an ellipsoid by a parallelotope
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
%    E = ellipsoid.generateRandom('Dimension',2);
%    Z = enc_parallelotope(E);
%
%    figure; hold on;
%    plot(E,[1,2],'b');
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

Tinv = sqrtm(E.Q);
n = length(E.Q);
% transform ellipsoid into sphere -> square around sphere -> back transform
Z = zonotope([E.q,Tinv*eye(n)]);

%------------- END OF CODE --------------