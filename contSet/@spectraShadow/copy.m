function SpS_out = copy(SpS)
% copy - copies the spectraShadow object (used for dynamic dispatch)
%
% Syntax:
%    SpS_out = copy(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    SpS_out - copied spectraShadow object
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    SpS_out = copy(SpS);
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
SpS_out = spectraShadow(SpS);

% ------------------------------ END OF CODE ------------------------------
