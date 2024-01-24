function cPZ = empty(n)
% empty - instantiates an empty constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    cPZ - empty constrained polynomial zonotope
%
% Example: 
%    cPZ = conPolyZono.empty(2);
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

cPZ = conPolyZono(zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
