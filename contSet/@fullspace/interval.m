function I = interval(fs)
% interval - converts a full-dimensional space to an interval object
%    case R^0: since the only contained point 0 is not representable in
%    MATLAB, we cannot convert R^0 to an interval
%
% Syntax:
%    I = interval(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    I - interval
%
% Example: 
%    fs = fullspace(2);
%    I = interval(fs);
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

if fs.dimension > 0
    % lower and upper bounds are infinity
    I = interval(-Inf(fs.dimension,1),Inf(fs.dimension,1));
else
    throw(CORAerror('CORA:outOfDomain','validDomain','>0'));
end

% ------------------------------ END OF CODE ------------------------------
