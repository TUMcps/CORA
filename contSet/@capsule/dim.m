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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       15-September-2019 
% Last update:   12-March-2021 (MW, add empty case)
%                09-January-2024 (MW, simplify)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = size(C.c,1);

% ------------------------------ END OF CODE ------------------------------
