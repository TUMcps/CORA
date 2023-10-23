function fs = project(fs,dims)
% project - projects a full-dimensional space onto the specified dimensions
%    case R^0: no dimensions for projection possible
%
% Syntax:
%    fs = project(fs,dims)
%
% Inputs:
%    fs - fullspace object
%    dims - dimensions for projection
%
% Outputs:
%    fs - projected fullspace
%
% Example: 
%    fs = fullspace(4);
%    val = project(fs,1:2);
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
    throw(CORAerror('CORA:notSupported','Projection of of R^0 not supported'));
elseif any(dims < 0 | dims > fs.dimension)
    throw(CORAerror('CORA:outOfDomain','validDomain',['1:' num2str(fs.dimension)]));
else
    fs.dimension = length(dims);
end

% ------------------------------ END OF CODE ------------------------------
