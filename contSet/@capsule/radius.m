function r = radius(C)
% radius - computes the radius of an enclosing hyperball of a capsule
%
% Syntax:
%    r = radius(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    r - radius of enclosing hyperball
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    r = radius(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-August-2019 
% Last update:   16-September-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute radius
r = norm(C.g,2) + C.r;

% ------------------------------ END OF CODE ------------------------------
