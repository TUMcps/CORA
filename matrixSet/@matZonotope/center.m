function c = center(matZ)
% center - Returns the center of an matZonotope
%
% Syntax:
%    c = center(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    c - center of the matrix zonotope
%
% Example:
%    M = matZonotope(eye(2),{eye(2),2*eye(2)});
%    c = center(M)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       23-July-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = matZ.center;

% ------------------------------ END OF CODE ------------------------------
