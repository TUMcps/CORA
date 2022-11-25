function plotPath(varargin);
% Purpose:  plot road object
% Pre:      road object
% Post:     ---
% Built:    04.12.06,MA
% Modified: 21.11.07,MA


%no color specified
if nargin==1
    obj=varargin{1}; %road object
    segments=[];
    plotStyle='k-'; %k=black

    
%segments defined
elseif nargin==2
    obj=varargin{1};
    segments=varargin{2};
    plotStyle='k-'; %k=black
    
%color defined
elseif nargin==3
    obj=varargin{1};
    segments=varargin{2};
    plotStyle=varargin{3};     
end

if isempty(segments)
    x=obj.segments.x;
    y=obj.segments.y;
else
    x=obj.segments.x(segments);
    y=obj.segments.y(segments);    
end
plot(x,y,plotStyle);