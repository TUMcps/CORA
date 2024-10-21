function fs_out = copy(fs)
% copy - copies the fullspace object (used for dynamic dispatch)
%
% Syntax:
%    fs_out = copy(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    fs_out - copied fullspace object
%
% Example: 
%    fs = fullspace(2);
%    fs_out = copy(fs);
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
fs_out = fullspace(fs);

% ------------------------------ END OF CODE ------------------------------
