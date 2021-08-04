function d = dim(C)
% dim - Returns the defined dimension of a capsule
%
% Syntax:  
%    d = dim(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    d - dimension of the capsule C
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    d = dim(C)
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
    d = length(center(C));
else
    d = 0;
end

%------------- END OF CODE --------------