function plot(matZ,varargin)
% plot - Plots 2-dimensional projection of a matrix zonotope
%
% Syntax:
%    plot(matZ)
%    plot(matZ,dimensions)
%    plot(matZ,dimensions,linespec)
%
% Inputs:
%    matZ - matZonotope object
%    dimensions - (optional) dimensions that should be projected
%    linespec - (optional) plot style
%
% Outputs:
%    -
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       22-June-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%convert from matrix zonotope to zonotope
Z=zonotope(matZ);
    
%plot zonotope
plot(Z,varargin{:});

% ------------------------------ END OF CODE ------------------------------
