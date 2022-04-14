function han = plot(probZ,varargin)
% plot - Plots 2-dimensional projection of a zonotope with a maximum
%    of 5 generators
%
% Syntax:  
%    h = plot(probZ) plots the probabilistic zonotope probZ for the first two dimensions
%    h = plot(probZ,dims) plots the probabilistic  zonotope probZ for the two dimensions i,j:
%        "dims=[i,j]" and returns handle to line-plot object
%    h = plot(probZ,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    probZ - probabilistic zonotope object
%    dims - (optional) dimensions that should be projected
%    type - (optional) plot settings (LineSpec and name-value pairs)
%    m - (optional) m-sigma value, default: probZ.gamma
%
% Outputs:
%    han - handle of graphics object
%    maxVal - maximum probabilistic value
%
% Example:
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    plot(pZ,'dark');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       03-August-2007
% Last update:   17-July-2020
% Last revision: ---

%------------- BEGIN CODE --------------

% default values
dims = [1,2];
filled = false; % coresponds to mesh/surf for prob. zonotopes
NVpairs = {};
m = probZ.gamma;

%If two arguments are passed    
if nargin==2
    dims=varargin{1};
    
%If three or more arguments are passed
elseif nargin>=3
    dims = varargin{1};
    % parse plot options
    [~,NVpairs] = readPlotOptions(varargin(2:end));
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
    [NVpairs,m] = readNameValuePair(NVpairs,'m','isscalar',m);
end 

%compute enclosing probability
eP = enclosingProbability(probZ,m,dims);

%plot and output the handle
if ~filled
    han = mesh(eP.X,eP.Y,eP.P,NVpairs{:});
else
    han = surf(eP.X,eP.Y,eP.P,NVpairs{:});
end
  

%------------- END OF CODE --------------