function Z = origin(n)
% origin - instantiates a zonotope that contains only the origin
%
% Syntax:
%    Z = zonotope.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    Z - zonotope representing the origin
%
% Example: 
%    Z = zonotope.origin(2);
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
Z = zonotope(zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
