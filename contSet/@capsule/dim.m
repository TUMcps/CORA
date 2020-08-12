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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

d = length(center(C));

%------------- END OF CODE --------------