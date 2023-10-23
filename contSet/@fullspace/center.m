function c = center(fs)
% center - returns the center of a full-dimensional space; we define this
%    to be the origin
%    case R^0: 0 (not representable in MATLAB)
%
% Syntax:
%    c = center(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    c - center
%
% Example: 
%    fs = fullspace(2);
%    c = center(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0
    c = NaN;
else
    c = zeros(fs.dimension,1);
end

% ------------------------------ END OF CODE ------------------------------
