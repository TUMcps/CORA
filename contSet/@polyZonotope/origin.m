function pZ = origin(n)
% origin - instantiates a polynomial zonotope that contains only the origin
%
% Syntax:
%    pZ = polyZonotope.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    pZ - polynomial zonotope representing the origin
%
% Example: 
%    pZ = polyZonotope.origin(2);
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
pZ = polyZonotope(zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
