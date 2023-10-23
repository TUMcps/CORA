function testRectangle(varargin)
% testRectangle - tests somne properties for the pre-computation of the 
%    fraction of a polytopes that intersects with another polytope of same 
%    size, but different position and orientation
%
% Syntax:
%    testRectangle()
%
% Inputs:
%    pos1 - ???
%    pos2 - ???
%    angle1 - ???
%    angle2 - ???
%
% Outputs:
%    fArray - 3-dimensional fraction array
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       02-April-2008 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin==0
    pos1=6*(rand(2,1)-0.5);
    pos2=6*(rand(2,1)-0.5);
    angle1=2*pi*(rand(1)-0.5);
    angle2=2*pi*(rand(1)-0.5);
else
    pos1=varargin{1};
    pos2=varargin{2};
    angle1=varargin{3};
    angle2=varargin{4};
end

%road data
segmentLength=6; %[m] %<-- enlarged due to curvature
roadWidth=4; %[m]
devSegments=7; 

%obtain length, width and including circle radius of the rectanular polytopes
length=segmentLength;
width=roadWidth/devSegments;
r=sqrt(length^2+width^2);

%compute beta angle (see sketch)
beta=atan(width/length);

%define deviation in the x-coordinate
deltaX=0.2; %in meters
deltaY=deltaX;

%determine optimal deltaAngle
deltaAngle=asin(deltaX/(cos(beta)*r));

%generate both polytopes; first polytope is not enlarged, second one is
Pol1=0.5*interval([-length;-width],[length;width]);
%enlarge first polytope
xEnlarge=deltaY+sin(beta)*sin(deltaAngle)*r;
yEnlarge=deltaX+cos(beta)*sin(deltaAngle)*r;
deltaPol=0.5*interval([-xEnlarge;-yEnlarge],[xEnlarge;yEnlarge]);
Pol2=Pol1+deltaPol;

%convert to zonotope
Z1=zonotope(Pol1);
Z2=zonotope(Pol2);

%generate rotation matrices
Rot1=[cos(angle1) -sin(angle1);...
     sin(angle1) cos(angle1)];
Rot2=[cos(angle2) -sin(angle2);...
     sin(angle2) cos(angle2)];

%generate 2 original zonotopes
Z1orig=Rot1*Z1+pos1;
Z2orig=Rot2*Z1+pos2;

%plot 2 original zonotopes
figure; 
hold on
axis equal
plot(Z1orig);
plot(Z2orig);


%generate rotation matrix
RotRel=[cos(-angle1) -sin(-angle1);...
     sin(-angle1) cos(-angle1)];

%return relAngle
relAngle=angle2-angle1;
relAngle=rem(relAngle,pi);

%return relPosition
relPos=RotRel*(pos2-pos1);

%generate positive position differences
if prod(relPos)<0
    relAngle=-relAngle;
end
relPos=abs(relPos);
%generate positive angle differences    
if relAngle<0
    relAngle=relAngle+pi;
end   
 
%generate rotation matrix
RotRel2=[cos(relAngle) -sin(relAngle);...
     sin(relAngle) cos(relAngle)];

%generate relative zonotope
Zrel=RotRel2*Z1+relPos;

%plot relative zonotopes
figure; 
hold on
axis equal
plot(Z1);
plot(Zrel);


%obtain angle and position segments
iAngleSeg=ceil(relAngle/deltaAngle);
iXseg=ceil(relPos(1)/deltaX);
iYseg=ceil(relPos(2)/deltaY);

%obtain x,y translations and angle rotations
xTrans=0.5*deltaX+(iXseg-1)*deltaX;
yTrans=0.5*deltaY+(iYseg-1)*deltaY;
angleTrans=0.5*deltaAngle+(iAngleSeg-1)*deltaAngle;
%generate rotation matrix
Rot=[cos(angleTrans) -sin(angleTrans);...
     sin(angleTrans) cos(angleTrans)];
 
%obtain enlarged zonotope
Zenl=Rot*Z2+[xTrans;yTrans];

%plot Zenl
plot(Zenl);
plot(Zrel);

%transform to polytopes and compute relative intersected volume
P1orig=polytope(Z1orig);
P2orig=polytope(Z2orig);
P1=polytope(Z1);
Prel=polytope(Zrel);
Penl=polytope(Zenl);

PintOrig=P1orig&P2orig;
PintRel=P1&Prel;
PintEncl=P1&Penl;

fullVol=modVolume(P1);

partVolOrig=modVolume(PintOrig)/fullVol
partVolRel=modVolume(PintRel)/fullVol
partVolEncl=modVolume(PintEncl)/fullVol

% ------------------------------ END OF CODE ------------------------------
