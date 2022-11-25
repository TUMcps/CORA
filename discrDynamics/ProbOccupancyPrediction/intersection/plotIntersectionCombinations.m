function plotIntersectionCombinations()
% plotIntersectionCombinations - plots the probability that two vehicle
% bodies intersect when their set of uncertain vehicle cemters is given
%
% Syntax:  
%    plotIntersectionCombinations()
%
% Inputs:
%    ---
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      22-October-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

figure
hold on

%obtain uncertain center interval
%road data
segmentLength=5; %[m] %<-- enlarged due to curvature
centerWidth=2; %[m]
devSegments=4; 

%obtain length, width and including circle radius of the rectanular polytopes
slength=segmentLength;
width=centerWidth/devSegments;

%obtain interval hull
IHcenter=0.5*interval([-slength; -width], [slength;width]);
IHcenterComb=0.5*interval([-slength; -centerWidth], [slength; centerWidth]);

%obtain uncertain body interval
%vehicle data
carLength=5; %[m] %<-- enlarged due to curvature
carWidth=2; %[m]

%obtain interval hull
IHcarOnly=0.5*interval([-carLength; -carWidth], [carLength; carWidth]);

%add both interval hulls
IHcar=IHcenter+IHcarOnly;
IHcarComb=IHcenterComb+IHcarOnly;

%obtain radius of interval hull
rCar=radius(IHcar);
rCarComb=radius(IHcarComb);

%compute polytope representation
Pcar=polytope(IHcarOnly);
P2=Pcar;

%set number of x, y and angle segments
nrOfxSeg=44;
nrOfySeg=44;
nrOfAngleSeg=72;


%define deviation in the x-coordinate
deltaX=rCar/nrOfxSeg; %in meters
deltaY=deltaX;
deltaAngle=pi/nrOfAngleSeg;

%set angle, x-segment, y-segment
iAngleSeg=50;
iXseg=23;
%iXseg=13;
iYseg=6;


%obtain x,y translations and angle rotations
xTrans=(iXseg-1)*deltaX;
yTrans=(iYseg-1)*deltaY;
angleTrans=(iAngleSeg-1)*deltaAngle;

%generate rotation matrix
Rot=[cos(angleTrans) -sin(angleTrans);...
     sin(angleTrans) cos(angleTrans)];

%modify second polytope
P2mod=Rot*P2+[xTrans;yTrans];

%compute grid of relative positions within the uncertain center
%regions
%relativePos=relativeGridPoints(IHcenter,Rot,[xTrans;yTrans],4);
relativePos=relativeGridPoints(IHcenter,Rot,[0;0],4);
nrOfRelativePos=max(size(relativePos(:,1)));

for iRelPos=1:nrOfRelativePos

    %intersect P and modified P2
    PcarInt=Pcar & (P2mod+relativePos(iRelPos,:)');
    
%     plot(Pcar);
%     plot(P2mod+relativePos(iRelPos,:)');

    %check if polytopes intersected
    %car
    [c,rad]=chebyball(PcarInt);
    if rad==-inf
        carInt(iRelPos)=0;
    else
        carInt(iRelPos)=1;
    end
end

%compute result
intersectionProbCar=max(size(find(carInt)))/nrOfRelativePos;


%convert interval hulls to zonotopes
Zfixed=zonotope(IHcenter);
Zvar=Rot*Zfixed+[xTrans;yTrans];

%plot uncertain position polytopes
plot(Zfixed,[1 2],'blackEdge');
plot(Zvar,[1 2],'blackEdge');


%plot grids of the polytopes
pointsFixed=gridPoints(IHcenter,4);
for i=1:length(pointsFixed(1,:))
    pointsVar(:,i)=Rot*pointsFixed(:,i)+[xTrans;yTrans];
end

%plot grid points
plot(pointsFixed(1,:),pointsFixed(2,:),'k+');
plot(pointsVar(1,:),pointsVar(2,:),'k+');

%plot exemplary vehicle bodies
extractedPointFixed=pointsFixed(:,8);
extractedPointVar=pointsVar(:,11);

%exemplary vehicle bodies
exemplaryBodyFixed=zonotope(IHcarOnly)+extractedPointFixed;
exemplaryBodyVar=Rot*zonotope(IHcarOnly)+extractedPointVar;

%plot exemplary vehicle bodies
plot(exemplaryBodyFixed,[1 2],'blackEdge');
plot(exemplaryBodyVar,[1 2],'blackEdge');

%plot centers
plot(extractedPointFixed(1),extractedPointFixed(2),'ko');
plot(extractedPointVar(1),extractedPointVar(2),'ko');

%------------- END OF CODE --------------