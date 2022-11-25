function [pos,angle]=segData(obj,segNumber,devNumber)
% segData - returns the position and angle of a segment for a given 
% segment number and deviation number.
%
% Syntax:  
%    [P]=segPolytope(obj,segNumber,devNumber)
%
% Inputs:
%    obj - road object
%    segNumber - number of the segment along the path
%    devNumber - number of the deviation segment
%
% Outputs:
%    P - segment polytope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 13-March-2008 
% Last update: 01-April-2008
% Last revision: ---

%------------- BEGIN CODE --------------


%get number of deviation segments
nrOfDev=obj.nrOfDevSegments;

%load data
angle=obj.segments.angle;
x=obj.segments.x;
y=obj.segments.y;


%obtain angles and positions at the beginning and end of the 
%path segment
%beginning of path segment:
angle1=angle(segNumber);
x1=x(segNumber);
y1=y(segNumber);
%end of path segment:
if segNumber==length(angle) %in case it is the last path segment
    angle2=angle1;
    x2=x1;
    y2=y1;    
else 
    angle2=angle(segNumber+1);
    x2=x(segNumber+1);
    y2=y(segNumber+1);    
end

%obtain middle point
xMid=x1+0.5*(x2-x1);
yMid=y1+0.5*(y2-y1);

%calculate translation vectors tangential to the path
transLat(1,1)=cos(angle1-0.5*pi)*obj.width/nrOfDev;
transLat(2,1)=sin(angle1-0.5*pi)*obj.width/nrOfDev;

%return position and angle
pos=[xMid; yMid]+(devNumber-0.5-0.5*nrOfDev)*transLat;
angle=angle1;

%------------- END OF CODE --------------