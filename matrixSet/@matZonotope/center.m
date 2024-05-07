function C = center(matZ)
% center - Returns the center of a matZonotope
%
% Syntax:
%    C = center(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    C - center of the matrix zonotope
%
% Example:
%    C = eye(2);
%    G(:,:,1) = eye(1);
%    G(:,:,2) = eye(2);
%    matZ = matZonotope(C,G);
%    C = center(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       23-July-2020 
% Last update:   ---
% Last revision: 25-April-2024 (TL, matZ.C property)

% ------------------------------ BEGIN CODE -------------------------------

C = matZ.C;

% ------------------------------ END OF CODE ------------------------------
