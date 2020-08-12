function points = pointSet(pZ,nrOfPoints)
% pointSet - Computes the set of random points inside the polynomial
%            zonotope
%
% Syntax:  
%    points = pointSet(pZ,nrOfPoints)
%
% Inputs:
%    pZ - polyZonotope object (n dimensional)
%    nrOfPoints - number of points in the resulting point cloud
%
% Outputs:
%    points - random point set (dimension: [n,M]) 
%
% Example: 
%    pZ = polyZonotope([0;0],[1 0;0 1],[0;0.1],[1,3]);
%
%    points = pointSet(pZ,100);
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(points(1,:),points(2,:),'.k'); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: pointSetExtreme

% Author:       Niklas Kochdumper
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% initialize point cloud
points = zeros(size(pZ.G,1),nrOfPoints);

% draw random points inside the polynomial zonotope
for i = 1:nrOfPoints
   points(:,i) = randPoint(pZ); 
end

%------------- END OF CODE --------------