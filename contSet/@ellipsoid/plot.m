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
%               12-March-2021
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

% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
end

% 2D vs 3D plot
if length(dims) == 2

    N = 1000;
    % project ellipsoid
    E_p = project(E,dims);
    %Since L only contains unit vectors l, we know that there exists a y
    %such that l'*y = suppfnc(E_p,L)=l'*q+sqrt(l'*Q*l). Therefore, y =
    %q+Q*l/sqrt(l'*Q*l).
%     if E_p.rank<=1
%         Y = boundary(E_p,2);%N=2 here returns exactly the two boundary points
%         % no filled in case of rank-1 ellipsoid
%         h = plot(Y(1,:),Y(2,:),linespec,NVpairs{:});
%         return;
%     end
    Y = boundary(E_p,N);
    %plot and output the handle
    if filled
        h = fill(Y(1,:),Y(2,:),linespec,NVpairs{:});
    else
        h = plot(Y(1,:),Y(2,:),linespec,NVpairs{:});
    end

else
    % enclose ellipsoid with zonotope
    E = project(E,dims);
    Z = zonotope(E,100,'o:norm');
    
    % plot zonotope enclosure
    h = plot(Z,[1,2,3],linespec,NVpairs{:},'Filled',filled);
end

%------------- END OF CODE --------------