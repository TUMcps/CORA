function res = isempty(pZ)
% isempty - checks if a polyZonotope object is empty
%
% Syntax:  
%    res = isempty(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - boolean whether pZ is empty
%
% Example: 
%    pZ1 = polyZonotope([2;1],[1 0; -2 1],[1; 0],[0 2; 1 0]);
%    isempty(pZ1);
%    pZ2 = polyZonotope([],[],[]);
%    isempty(pZ2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       01-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

res = isempty(pZ.c) && isempty(pZ.G) && isempty(pZ.Grest);

%------------- END OF CODE --------------

