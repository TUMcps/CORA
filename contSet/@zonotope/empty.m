function Z = empty(n)
% empty - instantiates an empty zonotope
%
% Syntax:
%    Z = zonotope.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    Z - empty zonotope
%
% Example: 
%    Z = zonotope.empty(2);
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

Z = zonotope(zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
