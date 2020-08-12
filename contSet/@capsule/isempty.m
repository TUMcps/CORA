function res = isempty(C)
% isempty - checks if capsule is empty
%
% Syntax:  
%    res = isempty(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    res - boolean whether capsule is empty or not
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    res = isempty(C) % false
%    C = capsule();
%    res = isempty(C); % true
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-Sep-2019
% Last update:  01-May-2020 (case r=0: empty)
% Last revision:---

%------------- BEGIN CODE --------------

res = isempty(center(C)) && isempty(C.g) && ...
    ( isempty(C.r) || (isscalar(C.r) && C.r == 0) );

%------------- END OF CODE --------------