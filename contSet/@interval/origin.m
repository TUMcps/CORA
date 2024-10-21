function I = origin(n)
% origin - instantiates an interval that contains only the origin
%
% Syntax:
%    I = interval.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    I - interval representing the origin
%
% Example: 
%    I = interval.origin(2);
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
I = interval(zeros(n,1),zeros(n,1));

% ------------------------------ END OF CODE ------------------------------
