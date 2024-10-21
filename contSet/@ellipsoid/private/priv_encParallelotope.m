function Z = priv_encParallelotope(E)
% priv_encParallelotope - encloses a non-degenerate ellipsoid by a
%    parallelotope
%
% Syntax:
%    Z = priv_encParallelotope(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[1;0]);
%    Z = zonotope(E,'outer:box',8);
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

% transform ellipsoid into sphere -> square around sphere -> back transform
Z = zonotope(E.q,sqrtm(E.Q));

% ------------------------------ END OF CODE ------------------------------
