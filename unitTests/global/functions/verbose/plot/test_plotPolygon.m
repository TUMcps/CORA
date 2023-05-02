function res = test_plotPolygon
% test_plotPolygon - unit test function for plotting of 2D/3D polygons
%
% Syntax:  
%    res = test_plotPolygon
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      30-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

try
    % 2D polygon
    V = [1 0; 1 2; 0 3; -2 2; -3 0; 0 -1; 1 0]';
    plotPolygon(V);
    plotPolygon(V,'r');
    plotPolygon(V,'FaceColor','r');

    % point cloud with convex hull
    V = randn(2,100);
    plotPolygon(V,'convHull',true);

    % 3D polygon
    V = [1 0 1; 1 2 -1; 0 3 0; -2 2 1; -3 0 2; 0 -1 0; 1 0 1]';
    plotPolygon(V);
    plotPolygon(V,'r');
    plotPolygon(V,'FaceColor','r');

    close all;

catch
    % close all figures
    close all;
    res = false;
end
%------------- END OF CODE --------------