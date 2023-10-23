function res = testSpecial_road_intersection
% testSpecial_road_intersection - unit_test_function for checking
% whether the intersection computation is correct
%
% Syntax:
%    res = testSpecial_road_intersection
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       18-August-2016
% Last update:   01-September-2016
%                12-July-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%load probabilistic model
load probModel_car_sim_02September2016_newPartitionFormat probModel
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
runs = 1e3;
runsMC = 1e3;

% init partial result
resPartial = [];

for i = 1:runs
    
    % display current run
    disp(num2str(i));

    %random generation of positions and orientations
    pos1 = randPoint(posZono);
    angle1 = randPoint(angleZono);
    pos2 = randPoint(posZono);
    angle2 = randPoint(angleZono);

    %obtain intersection probability from database (db)
    intersected_db = intersection_database(road(0,0,0),fArray,intDatabase,pos1,angle1,pos2,angle2);
    
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
    % intersection probbaility based on Monte Carlo simulation
    intersected_MC = sum(intArray)/runsMC;
    
    % init partial result
    resPartial(end+1) = 1;
    
    % there is a non-zero intersection probability although an intersection
    % is detected by the Monte carlo simulation
    if intersected_db == 0 && intersected_MC >= 0.1
        % test not passed
        resPartial(end) = 0;
        
        %plot for debugging
        aux_plotForDebugging(Zcenter, ZcenterEgo, bodyZono, ZcenterCurr, ...
                ZcenterEgoCurr, rect1, rect2, pos1, pos2, rot1, rot2);
    else
        relDiff = abs(intersected_db - intersected_MC)/max(intersected_MC,intersected_db);
        % the relative difference of the result from the database and the
        % Monte Carlo simulation should not exceed a threshold
        if intersected_db > 0.1 && intersected_MC > 0.1 && relDiff > 0.4
            % test not passed
            resPartial(end) = 0;
            
            %plot for debugging
            aux_plotForDebugging(Zcenter, ZcenterEgo, bodyZono, ZcenterCurr, ...
                ZcenterEgoCurr, rect1, rect2, pos1, pos2, rot1, rot2);
        end
    end
end

% compute final result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function aux_plotForDebugging(Zcenter, ZcenterEgo, bodyZono, ZcenterCurr, ...
    ZcenterEgoCurr, rect1, rect2, pos1, pos2, rot1, rot2)

    % plot situation
    figure 
    hold on
    plot(rot1*(Zcenter+bodyZono)+pos1,[1 2],'g:');
    plot(rot2*(ZcenterEgo+bodyZono)+pos2,[1 2],'g');
    plot(ZcenterCurr,[1 2],'r:');
    plot(ZcenterEgoCurr,[1 2],'r');
    draw(rect1,'b:');
    draw(rect2,'b');
end


% ------------------------------ END OF CODE ------------------------------
