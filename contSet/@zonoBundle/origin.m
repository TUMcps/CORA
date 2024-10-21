function zB = origin(n)
% origin - instantiates a zonotope bundle that contains only the origin
%
% Syntax:
%    zB = zonoBundle.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    zB - zonotope bundle representing the origin
%
% Example: 
%    zB = zonoBundle.origin(2);
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
zB = zonoBundle({zonotope(zeros(n,1))});

% ------------------------------ END OF CODE ------------------------------
