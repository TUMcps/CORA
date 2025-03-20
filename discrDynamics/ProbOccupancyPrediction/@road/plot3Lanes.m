function plot3Lanes(obj,except)
% plot3Lanes - plot road markings of AT08 evading maneuvers
%
% Syntax:
%    plot3Lanes(obj,except)
%
% Inputs:
%    road object
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
% Written:       24-April-2008
% Last update:   12-August-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

hold on

%load road variables
x=obj.segments.x;
y=obj.segments.y;
angle=obj.segments.angle;

for i=1:except(1)
    
        %obtain x and y coordinates
        x1=x(i);
        y1=y(i);
        angle1=angle(i);

        transLat(1)=cos(angle1-0.5*pi)*1*obj.width;
        transLat(2)=sin(angle1-0.5*pi)*1*obj.width;

        % points for left lane boundary
        xLeft(i)=x1-0.5*transLat(1);
        yLeft(i)=y1-0.5*transLat(2);

        % points for first lane divider
        xMid1(i)=x1+0.5*transLat(1);
        yMid1(i)=y1+0.5*transLat(2);    
        
        % points for second lane divider
        xMid2(i)=x1+1.5*transLat(1);
        yMid2(i)=y1+1.5*transLat(2);           

        % points for right lane boundary
        xRight(i)=x1+2.5*transLat(1);
        yRight(i)=y1+2.5*transLat(2);     
end

% plot individual lane boundaries
plot(xLeft,yLeft,'k-');
plot(xMid1,yMid1,'k--');
plot(xMid2,yMid2,'k--');
plot(xRight,yRight,'k-');

% ------------------------------ END OF CODE ------------------------------
