function han = plotRandPoint(pZ,varargin)
% plotRandPoint - Plots a point cloud of random points inside a polynomial
%                 zonotope
%
% Syntax:  
%    han = plot(pZ,dimensions,N,type)
%
% Inputs:
%    pZ - polyZonotope object
%    dimensions - dimensions that should be projected (optional) 
%    N - number of points for the point cloud (optional)
%    type - plot type (optional)
%
% Outputs:
%    han - handle for the resulting graphic object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%
%    hold on
%    plot(pZ,[1,2],'FaceColor',[.5 .5 .5],'Filled',true,'EdgeColor','none');
%    plotRandPoint(pZ,[1,2],100000,'.r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Author:       Nikals Kochdumper
% Written:      23-March-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dimensions = [1,2];
N = 1000;
type = {'.k'};

% parse input arguments
if nargin > 1
   dimensions = varargin{1}; 
end
if nargin > 2
   N = varargin{2}; 
end
if nargin > 3
   type = varargin(3:end);
end

% draw random points inside the polynomial zonotope
points = zeros(2,N);

for i = 1:N
   temp = randPoint(pZ);
   points(:,i) = [temp(dimensions(1));temp(dimensions(2))];
end

% plot the point cloud
han = plot(points(1,:),points(2,:),type{:});

%------------- END OF CODE --------------