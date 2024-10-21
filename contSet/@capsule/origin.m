function C = origin(n)
% origin - instantiates a capsule that contains only the origin
%
% Syntax:
%    C = capsule.origin(n)
%
% Inputs:
%    n - dimension (integer, >= 1)
%
% Outputs:
%    C - capsule representing the origin
%
% Example: 
%    C = capsule.origin(2);
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
C = capsule(zeros(n,1),zeros(n,1),0);

% ------------------------------ END OF CODE ------------------------------
