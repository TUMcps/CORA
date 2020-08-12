function [c] = center(pZ)
% center - Returns the center of a polynomial zonotope
%
% Syntax:  
%    [c] = center(qZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    c - center of the polynomial zonotope pZ
%
% Example: 
%    pZ = polyZonotope([3;4],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    center(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/center, quadZonotope/center

% Author:       Niklas Kochdumper
% Written:      27-June-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    c = pZ.c;

%------------- END OF CODE --------------