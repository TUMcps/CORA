function fs = Inf(n)
% Inf - instantiates a fullspace fullspace object
%
% Syntax:
%    fs = fullspace.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    fs - fullspace fullspace object
%
% Example: 
%    fs = fullspace.Inf(2);
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

% call constructor (and check n there)
fs = fullspace(n);

% ------------------------------ END OF CODE ------------------------------
