function pZ_out = copy(pZ)
% copy - copies the polyZonotope object (used for dynamic dispatch)
%
% Syntax:
%    pZ_out = copy(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    pZ_out - copied polyZonotope object
%
% Example: 
%    pZ = polyZonotope([2;1],[1 0; -2 1],[1;0],[0 2; 1 0]);
%    pZ_out = copy(pZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
pZ_out = polyZonotope(pZ);

% ------------------------------ END OF CODE ------------------------------
