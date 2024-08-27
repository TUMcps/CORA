function han = plotRandPoint(S,varargin)
% plotRandPoint - Plots a point cloud of random points from inside a set
%
% Syntax:
%    han = plotRandPoint(S)
%    han = plotRandPoint(S,dims,N,type)
%
% Inputs:
%    pZ - contSet object
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

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       11-December-2022
% Last update:   29-July-2024 (TL, added 1d and 3d plotting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values for the optional input arguments
[dims,N,type] = setDefaultValues({[1,2],1000,'.k'},varargin);

% check input arguments
inputArgsCheck({{S,'att','contSet'};
                {dims,'att','numeric',{'vector','integer','positive'}};
                {N,'att','numeric',{'positive','integer','scalar'}};
                {type,'att','char'}});

% generate N random points inside the set
points = randPoint_(S,N,'standard');

% plot the point cloud projected onto desired dimensions
if isscalar(dims)
    han = plot(points(dims(1),:),0*points(dims(1),:),type);
elseif numel(dims) == 2
    han = plot(points(dims(1),:),points(dims(2),:),type);
elseif numel(dims) == 3
    han = plot3(points(dims(1),:),points(dims(2),:),points(dims(3),:),type);
else
    throw(CORAerror('CORA:wrongValue','second','Number of dimensions to plot has to be <= 3.'))
end

% only return handle if desired
if nargout == 0
    clear han;
end

% ------------------------------ END OF CODE ------------------------------
