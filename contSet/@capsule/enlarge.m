function C = enlarge(C,factor)
% enlarge - enlarges the capsule around its center
%
% Syntax:  
%    C = enlarge(C,factor)
%
% Inputs:
%    C - capsule
%    factor - enlargement factor (scalar)
%
% Outputs:
%    C - capsule
%
% Example: 
%    C = capsule([1; 1], [0.5; -1], 0.5);
%    factor = 2;
%    C = enlarge(C,factor)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

C = capsule(center(C),factor*C.g,factor*C.r);

%------------- END OF CODE --------------