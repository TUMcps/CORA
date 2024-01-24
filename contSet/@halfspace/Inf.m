function hs = Inf(n)
% Inf - instantiates a fullspace halfspace object
%
% Syntax:
%    hs = halfspace.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    hs - fullspace halfspace object
%
% Example: 
%    hs = halfspace.Inf(2);
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

% halfspace 0*x <= 1 is fulfilled for all values for x
hs = halfspace(zeros(1,n),1);

% ------------------------------ END OF CODE ------------------------------
