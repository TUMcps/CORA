function n = dim(C)
% dim - returns the dimension of the ambient space of a capsule
%
% Syntax:  
%    n = dim(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    n - dimension of the ambient space
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    n = dim(C)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      15-Sep-2019 
% Last update:  12-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

if ~isempty(C)
    n = length(center(C));
else
    n = 0;
end

%------------- END OF CODE --------------