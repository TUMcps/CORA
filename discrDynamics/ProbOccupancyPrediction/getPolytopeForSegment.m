function [P] = getPolytopeForSegment(segments, segNumber, devNumber, nrOfDevSegments, widthOfLane)
%replaces segPolytope.m without using the CORA road class

%obtain angles and positions at the beginning and end of the path segment
%beginning of path segment:
angleStart = segments.angle(segNumber);
xStart = segments.x(segNumber);
yStart = segments.y(segNumber);

%end of path segment:
try %in case it is the last path segment
    angleEnd = segments.angle(segNumber+1);
    xEnd = segments.x(segNumber+1);
    yEnd = segments.y(segNumber+1);
catch
    angleEnd = angleStart;
    xEnd = xStart;
    yEnd = yStart;    
end

%calculate translation vectors tangential to the path
transLatStart(1,1) = cos(angleStart - (0.5 * pi)) * widthOfLane / nrOfDevSegments;
transLatStart(2,1) = sin(angleStart - (0.5 * pi)) * widthOfLane / nrOfDevSegments;

transLatEnd(1,1) = cos(angleEnd - (0.5 * pi)) * widthOfLane / nrOfDevSegments;
transLatEnd(2,1) = sin(angleEnd - (0.5 * pi)) * widthOfLane / nrOfDevSegments;

 
%define extreme points, l:left, r: right, b: bottom, t: top
Plb = [xStart; yStart] + (devNumber - 1 - (0.5 * nrOfDevSegments)) * transLatStart;
Prb = [xStart; yStart] + (devNumber - (0.5 * nrOfDevSegments)) * transLatStart;
Plt = [xEnd; yEnd] + (devNumber - 1 - (0.5 * nrOfDevSegments)) * transLatEnd;
Prt = [xEnd; yEnd] + (devNumber - (0.5 * nrOfDevSegments)) * transLatEnd;
V = [Plb, Prb, Plt, Prt];

%Generate polytope
P = polytope(V.');
end