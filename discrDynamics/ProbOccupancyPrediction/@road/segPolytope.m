function [P]=segPolytope(obj,segNumber,devNumber)
% segPolytope - returns the polytope of a reachable set segment where the
% segment number (along the path) and the deviation number (deviation from 
% the path) is given
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
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


%get number of deviation segments
nrOfDev=obj.nrOfDevSegments;


%obtain angles and positions at the beginning and end of the 
%path segment
%beginning of path segment:
angle1=obj.segments.angle(segNumber);
x1=obj.segments.x(segNumber);
y1=obj.segments.y(segNumber);
%end of path segment:
try %in case it is the last path segment
    angle2=obj.segments.angle(segNumber+1);
    x2=obj.segments.x(segNumber+1);
    y2=obj.segments.y(segNumber+1);
catch
    angle2=angle1;
    x2=x1;
    y2=y1;    
end

%calculate translation vectors tangential to the path
transLat1(1,1)=cos(angle1-0.5*pi)*obj.width/nrOfDev;
transLat1(2,1)=sin(angle1-0.5*pi)*obj.width/nrOfDev;

transLat2(1,1)=cos(angle2-0.5*pi)*obj.width/nrOfDev;
transLat2(2,1)=sin(angle2-0.5*pi)*obj.width/nrOfDev;

 
%define extreme points, l:left, r: right, b: bottom, t: top
Plb=[x1; y1]+(devNumber-1-0.5*nrOfDev)*transLat1;
Prb=[x1; y1]+(devNumber-0.5*nrOfDev)*transLat1;
Plt=[x2; y2]+(devNumber-1-0.5*nrOfDev)*transLat2;
Prt=[x2; y2]+(devNumber-0.5*nrOfDev)*transLat2;
V=[Plb,Prb,Plt,Prt];

%Generate polytope
%P=polytope(V.');
P=Polyhedron(V.');
    

%------------- END OF CODE --------------