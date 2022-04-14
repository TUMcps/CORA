function C = project(C,dims)
% project - Returns a capsule which is projected onto the specified
% dimensions
%
% Syntax:  
%    C = project(C,dims)
%
% Inputs:
%    C - capsule
%    dims - projected dimensions
%
% Outputs:
%    C - capsule
%
% Example: 
%    C = capsule([1; 1; 0], [0; 0; 1], 0.5);
%    Z = project(Z,[1 3]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%projection of center and generator; the radius is unchanged
C.c = C.c(dims);
if ~isempty(C.g)
    C.g = C.g(dims);
end

%------------- END OF CODE --------------