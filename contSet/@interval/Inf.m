function I = Inf(n)
% Inf - instantiates a fullspace interval
%
% Syntax:
%    I = interval.Inf(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    I - empty interval
%
% Example: 
%    I = interval.Inf(2);
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

I = interval(-Inf(n,1),Inf(n,1));

% ------------------------------ END OF CODE ------------------------------
