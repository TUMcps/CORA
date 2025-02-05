function plotPath(varargin)
% plotPath plot road object
% Pre:      road object
%
% Syntax:
%    plotPath(varargin)
%
% Inputs:
%    ???
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       04-December-2006
% Last update:   21-November-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
