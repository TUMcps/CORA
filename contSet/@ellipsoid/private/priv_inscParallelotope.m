function Z = priv_inscParallelotope(E)
% priv_inscParallelotope - inner-approximates a non-degenerate ellipsoid by
%    a parallelotope
%
% Syntax:
%    Z = priv_inscParallelotope(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,'inner:box',8);
% 
%    figure; hold on;
%    plot(E); plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       14-October-2019
% Last update:   27-January-2021 (MA, degenerate case implemented)
%                08-June-2021 (VG, moved degeneracy to main file)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(E);
T = inv(sqrtm(E.Q));
% transform ellipsoid into sphere -> square into sphere -> back transform
Z = zonotope(E.q,inv(T)*1/sqrt(n)*eye(n));

% ------------------------------ END OF CODE ------------------------------
