function P = polytope(fs)
% polytope - converts a full-dimensional space to a polytope object
%    case R^0: since the only contained point 0 is not representable in
%    MATLAB, we cannot convert R^0 to a polytope
%
% Syntax:
%    P = polytope(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    fs = fullspace(2);
%    P = polytope(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/Inf

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       14-December-2023
% Last update:   25-February-2025 (TL, used polytope.Inf)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if fs.dimension > 0
    % init polytope/inf with correct dimensions
    P = polytope.Inf(fs.dimension);
else
    throw(CORAerror('CORA:outOfDomain','validDomain','>0'));
end

% ------------------------------ END OF CODE ------------------------------
