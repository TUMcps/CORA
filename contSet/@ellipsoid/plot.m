function h=plot(E,varargin)
% plot - Plots 2-dimensional projection of an ellipsoid
%
% Syntax:  
%    h = plot(E) plots the ellipsoid E for the first two dimensions
%    h = plot(E,dims) plots the ellipsoid E for the two dimensions i,j:
%                   "dims=[i,j]" and returns handle to line-plot object
%    h = plot(E,dims,'Color','red',...) adds the standard plotting preferences
%
% Inputs:
%    E - ellipsoid object
%    dims - (optional) dimensions that should be projected
%    type - (optional) plot type (LineSpec and Name-Value pairs)
%
% Outputs:
%    handle
%
% Example: 
%    E=ellipsoid([1 0 0; 0 1 0;0 0 3]);
%    plot(E)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      13-March-2019
% Last update:  14-July-2020 (merge with plotFilled)
% Last revision:---

%------------- BEGIN CODE --------------

% default values
dims=[1,2];
linespec='b';
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

N = 1000;
t = linspace(0,2*pi,N);
L = [cos(t);sin(t)];
if length(E.Q)==1%scalar
    dims = 1;
end
% project ellipsoid
E_p = project(E,dims);
%Since L only contains unit vectors l, we know that there exists a y such
%that l'*y = suppfnc(E_p,L)=l'*q+sqrt(l'*Q*l). Therefore, y =
%q+Q*l/sqrt(l'*Q*l).
if length(E.Q)==1
    Y = boundary(E,2);%N=2 here returns exactly the two boundary points
    % no filled in case of one-dimensional ellipsoid
    h = plot(Y,zeros(size(Y)),linespec,NVpairs{:});
    set(gca,'ytick',[],'Ycolor','w','box','off')
    return;
end
Y = zeros(length(dims),size(L,2));
for i=1:size(L,2)
    l = L(:,i);
    Y(:,i) = E_p.q + E_p.Q*l/sqrt(l'*E_p.Q*l);
end
%plot and output the handle
if filled
    h = fill(Y(1,:),Y(2,:),linespec,NVpairs{:});
else
    h = plot(Y(1,:),Y(2,:),linespec,NVpairs{:});
end

%------------- END OF CODE --------------