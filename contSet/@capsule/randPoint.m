function p = randPoint(C)
% randPoint - Returns a random point of a capsule
%
% Syntax:  
%    p = randPoint(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    p - random point inside the capsule
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    p = randPoint(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

dimension = dim(C);
c = center(C);

randomDirection = rand(dimension,1);
randomDirection = randomDirection / norm(randomDirection,2);

p = c + (-1 + 2*rand(1)) * C.g + randomDirection * (rand(1)*C.r);


%------------- END OF CODE --------------