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
% Last update:   15-January-2024 (TL, parse input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    n = 0;
end
inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative'}}});

cZ = conZonotope(zeros(n,0),zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
