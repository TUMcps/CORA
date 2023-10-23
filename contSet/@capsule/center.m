function c = center(C)
% center - Returns the center of a capsule
%
% Syntax:
%    c = center(C)
%
% Inputs:
%    C - capsule
%
% Outputs:
%    c - center of the capsule C
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    c = center(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-March-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = C.c;

% ------------------------------ END OF CODE ------------------------------
