function cZ = origin(n)
% origin - instantiates a constrained zonotope that contains only the
%    origin
%
% Syntax:
%    cZ = conZonotope.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    cZ - constrained zonotope representing the origin
%
% Example: 
%    cZ = conZonotope.origin(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{n,'att','numeric',{'scalar','positive','integer'}}});
cZ = conZonotope(zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
