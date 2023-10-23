function fs = enclosePoints(points)
% enclosePoints - enclose a point cloud with a fullspace object
%
% Syntax:
%    fs = enclosePoints(points)
%
% Inputs:
%    points - point cloud (nxm matrix)
%
% Outputs:
%    fs - fullspace object
%
% Example: 
%    points = [2 4 -2; 1 0 5];
%    fs = fullspace.enclosePoints(points);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fs = fullspace(size(points,1));

% ------------------------------ END OF CODE ------------------------------
