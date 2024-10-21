function C_out = copy(C)
% copy - copies the capsule object (used for dynamic dispatch)
%
% Syntax:
%    C_out = copy(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    C_out - copied capsule object
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    C_out = copy(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
C_out = capsule(C);

% ------------------------------ END OF CODE ------------------------------
