function C = project(C,dims)
% project - projects a capsule onto the specified dimensions
%
% Syntax:
%    C = project(C,dims)
%
% Inputs:
%    C - (capsule) capsule
%    dims - dimensions for projection
%
% Outputs:
%    C - (capsule) projected capsule
%
% Example: 
%    C = capsule([1; 1; 0], [0; 0; 1], 0.5);
%    C1 = project(C,[1 3]);
%    C2 = project(C,[true false true]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%projection of center and generator; the radius is unchanged
C.c = C.c(dims);
if ~isempty(C.g)
    C.g = C.g(dims);
end

% ------------------------------ END OF CODE ------------------------------
