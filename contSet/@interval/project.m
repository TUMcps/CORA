function I = project(I,dims)
% project - projects an interval onto the specified dimensions
%
% Syntax:
%    res = project(I,dims)
%
% Inputs:
%    I - (interval) interval
%    dims - dimensions for projection
%
% Outputs:
%    I - (interval) projected interval
%
% Example:
%    I = interval([-3;-5;-2],[3;2;1]);
%    dims = [1,3];
%    project(I,dims)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   21-May-2022 (remove new instantiation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if all(size(I) > 1)
    % not implemented for matrices
    throw(CORAerror("CORA:wrongValue",'first','non-matrix interval object'));
else
    I = interval(I.inf(dims),I.sup(dims));
end

% ------------------------------ END OF CODE ------------------------------
