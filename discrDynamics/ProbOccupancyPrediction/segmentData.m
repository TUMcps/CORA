function [pos, angle] = segmentData(segments, segNumber, devNumber, nrOfDevSegments, widthOfLane)
% segmentData - returns the position and angle of a segment for a given 
% segment number and deviation number.
%
% Inputs:
%    segments - segments struct of vehicle path
%    segNumber - number of the segment along the path
%    devNumber - number of the deviation segment
%    nrOfDevSegments - the number of deviation segments
%
% Outputs:
%    P - segment polytope
%
% replaces segData.m without using the CORA road class


%load data
angle = segments.angle;
x = segments.x;
y = segments.y;

%obtain angles and positions at the beginning and end of the path segment
%beginning of path segment:
angle1 = angle(segNumber);
x1 = x(segNumber);
y1 = y(segNumber);
%end of path segment:
if segNumber == length(angle) %in case it is the last path segment
    angle2 = angle1;
    x2 = x1;
    y2 = y1;    
else 
    angle2 = angle(segNumber+1);
    x2 = x(segNumber+1);
    y2 = y(segNumber+1);    
end

%obtain middle point
xMid = x1 + 0.5 * (x2 - x1);
yMid = y1 + 0.5 * (y2 - y1);

%calculate translation vectors tangential to the path
transLat(1,1) = cos(angle1 - (0.5 * pi)) * widthOfLane / nrOfDevSegments;
transLat(2,1) = sin(angle1 - (0.5 * pi)) * widthOfLane / nrOfDevSegments;

%return position and angle
pos = [xMid; yMid] + (devNumber - 0.5 - (0.5 * nrOfDevSegments)) * transLat;
angle = angle1;

%------------- END OF CODE --------------