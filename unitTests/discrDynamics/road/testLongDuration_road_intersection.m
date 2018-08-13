function res = testLongDuration_road_intersection
%  testLongDuration_road_intersection - unit_test_function for checking
%  whether the intersection computation is correct
%
% Syntax:  
%    res = testLongDuration_road_intersection
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  01-September-2016
% Last revision:---


%------------- BEGIN CODE --------------

%set path
global filePath
filePath = [coraroot '/contDynamics/stateSpaceModels'];

%instantiate dummy road object
R = road(0,0,0);

%load probabilistic model
[FileName,PathName] = uigetfile();
cd(PathName);
file=load(FileName);
probModel=file.probModel;
fArray = probModel.fArray;

%extract required information
intDatabase = fArray.val.car;
rect = fArray.R.car;
%store road information
segLength = fArray.segLengthOther;
segWidth = fArray.roadWidth/fArray.devSegments;
segLengthEgo = fArray.segLengthEgo;
segWidthEgo = fArray.segWidthEgo;

%obtain interval of center uncertainty of other traffic participants
Icenter = 0.5*interval([-segLength;-segWidth], [segLength;segWidth]);
Zcenter = zonotope(Icenter);

%obtain interval of center uncertainty of ego vehicle
IcenterEgo = 0.5*interval([-segLengthEgo;-segWidthEgo], [segLengthEgo;segWidthEgo]);
ZcenterEgo = zonotope(IcenterEgo);

posInterval = interval([-15;-15],[15;15]);
posZono = zonotope(posInterval);
angleInterval = interval(-pi,pi);
angleZono = zonotope(angleInterval);
bodyInterval = 0.5*interval([-rect.width;-rect.height],[rect.width;rect.height]);
bodyZono = zonotope(bodyInterval);

%set number of runs
runs = 1e4;
runsMC = 100;

for i = 1:runs

    %random generation of positions and orientations
    pos1 = randPoint(posZono);
    angle1 = randPoint(angleZono);
    pos2 = randPoint(posZono);
    angle2 = randPoint(angleZono);

    %obtain intersection probability from database (db)
    intersected_db = intersection_database(R,fArray,intDatabase,pos1,angle1,pos2,angle2);
    
    %obtain rotation matrices
    rot1 = [cos(angle1) -sin(angle1);...
            sin(angle1) cos(angle1)];
    rot2 = [cos(angle2) -sin(angle2);...
            sin(angle2) cos(angle2)];
        
    %obtain uncertain centers
    ZcenterCurr = rot1*Zcenter + pos1;
    ZcenterEgoCurr = rot2*ZcenterEgo + pos2;
    
    %compute intersection probability based on Monte Carlo
    for iMC=1:runsMC
        % obtain random center points
        randCenterPos = randPoint(ZcenterCurr);
        randCenterPosEgo = randPoint(ZcenterEgoCurr);
        
        % obtain current rectangles
        rect1 = Rectangle(rect.width, rect.height, angle1, randCenterPos);
        rect2 = Rectangle(rect.width, rect.height, angle2, randCenterPosEgo);
        
        % check if rectangles intersect
        intArray(iMC) = intersect(rect1, rect2);
    end
    
    intersected_MC = sum(intArray)/runsMC;
    
    if intersected_db == 0
        if intersected_MC >= 0.05
            figure 
            hold on
            plot(rot1*(Zcenter+bodyZono)+pos1,[1 2],'g:');
            plot(rot2*(ZcenterEgo+bodyZono)+pos2,[1 2],'g');
            plot(ZcenterCurr,[1 2],'r:');
            plot(ZcenterEgoCurr,[1 2],'r');
            draw(rect1,'b:');
            draw(rect2,'b');
            
            %plot equivalent situation
            figure 
            hold on
            rot1neg = [cos(-angle1) -sin(-angle1);...
            sin(-angle1) cos(-angle1)];
            relAngle=angle2-angle1;
            rotDiff = [cos(relAngle) -sin(relAngle);...
            sin(relAngle) cos(relAngle)];
            relPos=rot1neg*(pos2-pos1);
            plot((Zcenter+bodyZono),[1 2],'g:');
            plot(rotDiff*(ZcenterEgo+bodyZono)+relPos,[1 2],'g');
            
            %change to only positive values
            relAngle=rem(relAngle,pi);
            if prod(relPos)<0
                relAngle=-relAngle;
            end
            relPos=abs(relPos);
            if relAngle<0
                relAngle=relAngle+pi;
            end 
            rotDiff = [cos(relAngle) -sin(relAngle);...
            sin(relAngle) cos(relAngle)];
            plot(rotDiff*(ZcenterEgo+bodyZono)+relPos,[1 2],'b');
        end
    else
        relDiff = abs(intersected_db - intersected_MC)/max(intersected_MC,intersected_db);
        if relDiff > 0.2
            figure 
            hold on
            plot(rot1*(Zcenter+bodyZono)+pos1,[1 2],'g:');
            plot(rot2*(ZcenterEgo+bodyZono)+pos2,[1 2],'g');
            plot(ZcenterCurr,[1 2],'r:');
            plot(ZcenterEgoCurr,[1 2],'r');
            draw(rect1,'b:');
            draw(rect2,'b');
            
            %plot equivalent situation
            figure 
            hold on
            rot1neg = [cos(-angle1) -sin(-angle1);...
            sin(-angle1) cos(-angle1)];
            relAngle=angle2-angle1;
            rotDiff = [cos(relAngle) -sin(relAngle);...
            sin(relAngle) cos(relAngle)];
            relPos=rot1neg*(pos2-pos1);
            plot((Zcenter+bodyZono),[1 2],'g:');
            plot(rotDiff*(ZcenterEgo+bodyZono)+relPos,[1 2],'g');
            
            %change to only positive values
            relAngle=rem(relAngle,pi);
            if prod(relPos)<0
                relAngle=-relAngle;
            end
            relPos=abs(relPos);
            if relAngle<0
                relAngle=relAngle+pi;
            end 
            rotDiff = [cos(relAngle) -sin(relAngle);...
            sin(relAngle) cos(relAngle)];
            plot(rotDiff*(ZcenterEgo+bodyZono)+relPos,[1 2],'b');
        end
    end
end





%------------- END OF CODE --------------
