function O_out = copy(O)
% copy - copies the emptySet object (used for dynamic dispatch)
%
% Syntax:
%    O_out = copy(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    O_out - copied emptySet object
%
% Example: 
%    O = emptySet(2);
%    O_out = copy(O);
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
O_out = emptySet(O);

% ------------------------------ END OF CODE ------------------------------
