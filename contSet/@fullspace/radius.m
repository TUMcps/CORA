function val = radius(fs)
% radius - computes the radius of the enclosing hyperball of a
%    full-dimensional space
%    case R^0: NaN
%
% Syntax:
%    val = radius(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    val - radius
%
% Example: 
%    fs = fullspace(2);
%    val = radius(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   25-April-2023 (MW, integrate R^0 case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension == 0
    val = NaN;
else
    val = Inf;
end

% ------------------------------ END OF CODE ------------------------------
