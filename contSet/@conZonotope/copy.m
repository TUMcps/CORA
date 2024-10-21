function cZ_out = copy(cZ)
% copy - copies the conZonotope object (used for dynamic dispatch)
%
% Syntax:
%    cZ_out = copy(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    cZ_out - copied conZonotope object
%
% Example: 
%    c = [0;0]; G = [1 0 1;0 1 1];
%    A = [1 1 1]; b = 4;
%    cZ = conZonotope([c,G],A,b);
%    cZ_out = copy(cZ);
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
cZ_out = conZonotope(cZ);

% ------------------------------ END OF CODE ------------------------------
