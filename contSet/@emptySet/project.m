function O = project(O,dims)
% project - projects an empty set onto the specified dimensions
%
% Syntax:
%    O = project(O,dims)
%
% Inputs:
%    O - emptySet object
%    dims - dimensions for projection
%
% Outputs:
%    O - projected emptySet object
%
% Example: 
%    O = emptySet(4);
%    project(O,[1,3])
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

if any(dims < 0 | dims > O.dimension)
    throw(CORAerror('CORA:outOfDomain','validDomain',['1:' num2str(O.dimension)]));
else
    O.dimension = length(dims);
end

% ------------------------------ END OF CODE ------------------------------
