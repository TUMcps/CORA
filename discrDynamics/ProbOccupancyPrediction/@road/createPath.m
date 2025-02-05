function [Obj]=createPath(varargin)
% createPath - ???
%
% Syntax:
%    [Obj]=createPath(varargin)
%
% Inputs:
%    ???
%
% Outputs:
%    ???
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       ???
% Written:       ???
% Last update:   14-November-2007
%                29-August-2008
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%standard specification
if nargin==4
    Obj=varargin{1};
    init=varargin{2};
    deltaAngle=varargin{3};
    timeSteps=varargin{4};
    sL=Obj.segmentLength;

%segment length defined
elseif nargin==5
    Obj=varargin{1};
    init=varargin{2};
    deltaAngle=varargin{3};
    timeSteps=varargin{4};
    sL=varargin{5};
end

%init
angle(1)=init(1);
x(1)=init(2);
y(1)=init(3);
j=1;

for i=1:(length(deltaAngle))
    for iStep=1:timeSteps(i)
        angle(j+1)=angle(j)+deltaAngle(i);
        %get x-values
        x(j+1)=x(j)+sL*cos(angle(j));
        %get y-values
        y(j+1)=y(j)+sL*sin(angle(j));    
        %counter
        j=j+1;
    end
end
%plot(x,y);
%write to object
Obj.segments.x=x;
Obj.segments.y=y;
Obj.segments.angle=angle;

save('path.mat', 'x', 'y', 'angle');

% ------------------------------ END OF CODE ------------------------------
