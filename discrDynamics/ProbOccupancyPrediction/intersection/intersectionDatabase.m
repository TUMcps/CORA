function fArray =intersectionDatabase(parameters)
% intersectionDatabase - pre-computes the fraction of a polytopes that
% intersects with another polytope of same size, but different position and
% and orientation
%
% Syntax:  
%    P = intersectionDatabase(parameters)
%
% Inputs:
%    parameters - struct of parameters mainly regarding cell sizes
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

% Author:       Matthias Althoff
% Written:      01-April-2008 
% Last update:  28-January-2009
%               17-June-2009
%               22-October-2009
%               01-August-2016
%               03-August-2016
%               18-August-2016
%               01-August-2017
% Last revision: ---

%------------- BEGIN CODE --------------


%obtain uncertain center interval
%other traffic participants
segLength = parameters.segLength; %[m] 
centerWidth = parameters.centerWidth; %[m]
devSegments = parameters.devSegments; 
segWidth = centerWidth/devSegments;

%ego vehicle
segLengthEgo = parameters.segLengthEgo; %[m] 
segWidthEgo = parameters.segWidthEgo; %[m]

%obtain uncertain body interval
%vehicle data
carLength = parameters.carLength; %[m] 
carWidth = parameters.carWidth; %[m]

bicycleLength = parameters.bicycleLength; %[m] 
bicycleWidth = parameters.bicycleWidth; %[m]

%set number of x, y and angle segments
nrOfxSeg = parameters.nrOfxSeg;
nrOfySeg = parameters.nrOfySeg;
nrOfAngleSeg = parameters.nrOfAngleSeg;

%obtain interval of center uncertainty of other traffic participants
IHcenter = 0.5*interval([-segLength;-segWidth], [segLength;segWidth]);

%obtain interval of center uncertainty of ego vehicle
IHcenterEgo = 0.5*interval([-segLengthEgo;-segWidthEgo], [segLengthEgo;segWidthEgo]);

%enclosing rectangle of car and bicycle
Rcar = Rectangle(carLength,carWidth,0,[0, 0]);
Rbicycle = Rectangle(bicycleLength,bicycleWidth,0,[0, 0]);

%enclosing rectangle of car and bicycle plus uncertain center
Rcar_unc = Rectangle(carLength+segLength,carWidth+segWidth,0,[0, 0]);
Rbicycle_unc = Rectangle(bicycleLength+segLength,bicycleWidth+segWidth,0,[0, 0]);

%enclosing rectangle of car and bicycle plus combined uncertain center
Rcar_uncComb = Rectangle(carLength+segLength,carWidth+centerWidth,0,[0, 0]);
Rbicycle_uncComb = Rectangle(bicycleLength+segLength,bicycleWidth+centerWidth,0,[0, 0]);

%enclosing rectangles of ego vehicle
Rego_unc = Rectangle(carLength+segLengthEgo,carWidth+segWidthEgo,0,[0, 0]);

%obtain radius of interval hull
rCar = getradius(Rcar_unc);
rBicycle = getradius(Rbicycle_unc);
rEgo = getradius(Rego_unc);
rCarComb = getradius(Rcar_uncComb);
rBicycleComb = getradius(Rbicycle_uncComb);


%define deviation in the x-coordinate
deltaX = (rEgo + rCarComb)/nrOfxSeg; %in meters
deltaY = deltaX;
deltaAngle = pi/nrOfAngleSeg;

for iAngleSeg = 1:nrOfAngleSeg
    for iXseg = 1:nrOfxSeg
        iAngleSeg
        iXseg
        for iYseg = 1:nrOfySeg
            %obtain x,y translations and angle rotations
            xTrans = (iXseg-1)*deltaX;
            yTrans = (iYseg-1)*deltaY;
            angleTrans = (iAngleSeg-1)*deltaAngle;
            
            %compute grid of relative positions within the uncertain center
            %regions
            %other vehicle
            relativePos = relativeGridPoints(IHcenter,0,[0;0],4);
            nrOfRelativePos = length(relativePos(1,:));
            
            %ego vehicle
            relativePosEgo = relativeGridPoints(IHcenterEgo,angleTrans,[0;0],4);
            nrOfRelativePosEgo = length(relativePos(1,:));
            
%             figure
%             hold on
            
            %init
            carInt = [];
            bicycleInt = [];
            
            for iRelPosEgo=1:nrOfRelativePosEgo
                for iRelPos=1:nrOfRelativePos
                    %create enclosing rectangles
                    Rego_curr = Rectangle(carLength,carWidth,angleTrans,[xTrans + relativePosEgo(1,iRelPosEgo), yTrans + relativePosEgo(2,iRelPosEgo)]);
                    Rcar_curr = Rectangle(carLength,carWidth,0,[relativePos(1,iRelPos), relativePos(2,iRelPos)]);
                    Rbicycle_curr = Rectangle(bicycleLength,bicycleWidth,0,[relativePos(1,iRelPos), relativePos(2,iRelPos)]);
                    
%                     figure
%                     hold on
%                     draw(Rcar_curr,'g');
%                     draw(Rego_curr,'r');

                    %check if rectangles intersect
                    carInt(end+1) = intersect(Rego_curr, Rcar_curr);
                    bicycleInt(end+1) = intersect(Rego_curr, Rbicycle_curr);
                end
            end
            
            %save result
            intersectionCar(iXseg,iYseg) = length(find(carInt))/(nrOfRelativePos*nrOfRelativePosEgo);
            intersectionBicycle(iXseg,iYseg) = length(find(bicycleInt))/(nrOfRelativePos*nrOfRelativePosEgo);
        end
    end
    %save results
    fArray.val.car{iAngleSeg} = sparse(intersectionCar);
    fArray.val.bicycle{iAngleSeg} = sparse(intersectionBicycle); 
end

%store road information
fArray.segLengthOther = segLength;
fArray.roadWidth = centerWidth; 
fArray.devSegments = devSegments;
fArray.segLengthEgo = segLengthEgo;
fArray.segWidthEgo = segWidthEgo;
%store segment lengths
fArray.segLength.angle = deltaAngle;
fArray.segLength.x = deltaX;
fArray.segLength.y = deltaY;
%store radii
fArray.radius.r.car = rCar;
fArray.radius.r.bicycle = rBicycle;
fArray.radius.rComb.car = rCarComb;
fArray.radius.rComb.bicycle = rBicycleComb;
fArray.radius.r.ego = rEgo;

%store rectangles
fArray.R.car = Rcar;
fArray.R.bicycle = Rbicycle;

%------------- END OF CODE --------------


function [movedPoints]=relativeGridPoints(obj,angleTrans,trans,intervals)
% relativeGridPoints - Computes grid points of two identical interval 
% hulls using gridPoints(); the second hull is translated and rotated and
% all combinations of relative positions of grid points are returned
%
% Syntax:  
%    [coordinateMat]=relativeGridPoints(obj,angleTrans,trans,intervals)
%
% Inputs:
%    obj - interval hull object
%    rot - rotation matrix
%    trans - translation vector
%    intervals - number of intervals for each dimension
%
% Outputs:
%    coordinateMat - matrix where columns are the relative positions of
%    grid points
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      17-June-2009 
% Last update:  03-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute grid points
% originalPoints=gridPoints(obj,intervals);
originalPoints=aux_gridPoints(obj,intervals);

%generate rotation matrix
rot=[cos(angleTrans) -sin(angleTrans);...
     sin(angleTrans) cos(angleTrans)];

%compute translated and rotated points
nrOfPoints=length(originalPoints(1,:));
movedPoints=rot*originalPoints+trans*ones(1,nrOfPoints);


function coordinateMat = aux_gridPoints(obj,segments)

%obtain segment length
segLengthVec = (obj.sup-obj.inf)./segments;

%obtain first grid point
startingPoint = obj.inf(:,1) + 0.5*segLengthVec;
coordinateMat = startingPoint;

%obtain segment combinations
n = length(obj.inf(:,1));

if isscalar(segments)
    comb = full_fact_mod(ones(n,1)*segments);
else
    comb = full_fact_mod(segments);
end

%add further grid points
for i = 2:length(comb(:,1))
    coordinateMat(:,i) = startingPoint+(comb(i,:)'-1).*segLengthVec;
end

%------------- END OF CODE --------------