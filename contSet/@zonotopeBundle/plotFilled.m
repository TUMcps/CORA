function handle = plotFilled(varargin)
% plotFilled - Plots 2-dimensional over-approximative projection of a 
% zonotope bundle and fills it out
%
% Syntax:  
%    handle = plot(Z,dimensions,plotOpt)
%
% Inputs:
%    Z - zonotope bundle object
%    dimensions - dimensions that should be projected (optional)
%    plotOpt - plot options (color, edge color, etc.)
%
% Outputs:
%    handle - handle to the graphics object
%
% Example: 
%    Z{1} = zonotope([0 1 1 -1;0 0 2 1]);
%    Z{2} = zonotope([2 1 1 -1;1 0 -1 -1]);
%    zB = zonotopeBundle(Z);
%
%    plotFilled(zb,[1,2],'r','EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      18-August-2016
% Last update:  19-October-2015 (NK, accelerate by using polyshape class)
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin==1
    Zbundle=varargin{1};
    dimensions=[1,2];
    plotOpt={'b'};
    
%If two arguments are passed    
elseif nargin==2
    Zbundle=varargin{1};
    dimensions=varargin{2};
    plotOpt={'b'};
   
%If three or more arguments are passed
elseif nargin>=3
    Zbundle=varargin{1};
    dimensions=varargin{2};   
    plotOpt=varargin(3:end);
end

w = warning;
warning('off','all');

%compute polytopes
for i=1:Zbundle.parallelSets
    %delete zero generators
    Z=deleteZeros(Zbundle.Z{i});
    %project zonotope
    Z=project(Z,dimensions);
    %convert to polytope
    temp = polygon(Z);
    V = temp(:,2:end);
    P{i} = polyshape(V(1,:),V(2,:));
end

%intersect polytopes
Pint=P{1};
for i=2:Zbundle.parallelSets
    Pint = intersect(Pint,P{i});
end

warning(w);
    
%Plot polytope
handle = plot(Pint,'FaceAlpha',1,'FaceColor',plotOpt{:});

%------------- END OF CODE --------------