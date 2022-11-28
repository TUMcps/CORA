function han = plotRandPoint(pZ,varargin)
% plotRandPoint - Plots a point cloud of random points from inside a
%    polynomial zonotope
%
% Syntax:  
%    han = plot(pZ,dims,N,type)
%
% Inputs:
%    pZ - polyZonotope object
%    dims - (optional) dimensions that should be projected
%    N - (optional) number of points for the point cloud
%    type - (optional) plot type
%
% Outputs:
%    han - handle for the resulting graphics object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor',[.5 .5 .5]);
%    plotRandPoint(pZ);
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
[dims,N,type] = setDefaultValues({[1,2],1000,{'.k'}},varargin{:});

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {dims,'att','numeric','nonnan'};
                {N,'att','numeric','nonnegative'};
                {type,'att','cell','nonempty'}});

% draw random points inside the polynomial zonotope
points = zeros(2,N);

for i = 1:N
    temp = randPoint(pZ);
    points(:,i) = [temp(dims(1));temp(dims(2))];
end

% plot the point cloud
han = plot(points(1,:),points(2,:),type{:});

if nargout == 0
    clear han;
end

%------------- END OF CODE --------------