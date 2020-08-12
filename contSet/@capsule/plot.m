function han = plot(C,varargin)
% plot - Plots 2-dimensional projection of a capsule
%
% Syntax:  
%    h = plot(C) plots the capsule for the first two dimensions
%    h = plot(C,dims) plots the capsule for the two dimensions i,j:
%                   "dims=[i,j]" and returns handle to line-plot object
%    h = plot(C,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    C - capsule object
%    dims - (optional) dimensions that should be projected
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to graphics object
%
% Example: 
%    C = capsule([1; 1; 0], [0; 1; 1], 0.5);
%    plot(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default
dims=[1,2];
linespec = 'b';
filled = false;
NVpairs = {};
    
%If two arguments are passed    
if nargin==2
    dims=varargin{1};
    
%If three or more arguments are passed
elseif nargin>=3
    dims=varargin{1};   
    % parse plot options
    [linespec,NVpairs] = readPlotOptions(varargin(2:end));
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
end

% project zonotope
C = project(C,dims);

% underapproximate capsule by polygon
p = polygon(C);

%plot and output the handle
if filled
    han = fill(p(1,:),p(2,:),linespec,NVpairs{:});
else
    han = plot(p(1,:),p(2,:),linespec,NVpairs{:});
end

%------------- END OF CODE --------------