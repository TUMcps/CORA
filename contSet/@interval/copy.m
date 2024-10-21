function I_out = copy(I)
% copy - copies the interval object (used for dynamic dispatch)
%
% Syntax:
%    I_out = copy(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I_out - copied interval object
%
% Example: 
%    I = interval([-2;1],[1;2]);
%    I_out = copy(I);
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
I_out = interval(I);

% ------------------------------ END OF CODE ------------------------------
