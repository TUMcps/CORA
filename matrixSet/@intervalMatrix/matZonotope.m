function matZ = matZonotope(intMat)
% matZonotope - converts an interval matrix to a matrix zonotope
%
% Syntax:
%    matZ = matZonotope(intMat)
%
% Inputs:
%    intMat - intervalMatrix object
%
% Outputs:
%    matZ - matrix zonotope
%
% Example: 
%    c = [2 3 4; 5 6 7];
%    d = [1 0 1; 0 0 1];
%    intMat = intervalMatrix(c,d);
%    matZ = matZonotope(intMat);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger, Tobias Ladner
% Written:       21-June-2010 
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                18-June-2023 (MW, fix conversion of interval matrices(!))
%                25-April-2024 (TL, much faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to interval
I = interval(intMat);

% read out dimensions
dims = dim(intMat);
n = dims(1); m = dims(2);
h = n*m;

% read out center
C = center(I);

% compute radius
r = rad(I);

% get generators
idx = 1:(n*m+1):n*m*h;
G = zeros(n,m,h);
G(idx) = r;

% instantiate matrix zonotope
matZ = matZonotope(C,G);

% ------------------------------ END OF CODE ------------------------------
