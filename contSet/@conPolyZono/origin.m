function cPZ = origin(n)
% origin - instantiates a constrained polynomial zonotope that contains
%    only the origin
%
% Syntax:
%    cPZ = conPolyZono.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    cPZ - conPolyZono representing the origin
%
% Example: 
%    cPZ = conPolyZono.origin(2);
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
cPZ = conPolyZono(zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
