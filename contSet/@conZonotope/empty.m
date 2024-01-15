function cZ = empty(n)
% empty - instantiates an empty constrained zonotope
%
% Syntax:
%    cZ = conZonotope.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    cZ - empty constrained zonotope
%
% Example: 
%    cZ = conZonotope.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

cZ = conZonotope(zeros(n,0),zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
