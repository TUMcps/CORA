function [statObs,dynObs,x0,goalSet,lanelets, information] = commonroad2cora(filename, varargin)
% commonroad2cora - convert CommonRoad XML-file scenario to CORA model
%
% Syntax:
%    [statObs,dynObs,x0,goalSet,lanelets] = commonroad2cora(filename);
%
% Inputs:
%    filename - name of the CommonRoad XML-file scenario (specified as string)
%    'verbose', true/false - name/input pair specifying whether conversion
%           information should be printed to the terminal. Default is true.
%
% Outputs:
%    statObs - cell-array storing the static obstacles
%    dynObs - cell-array storing the set and the time interval for the
%             dynamic obstacles
%    x0 - initial point for the planning problem
%    goalSet - cell-array storing the goal sets and the corresponding time
%    lanelets - cell-array storing the sets for all lanelets
%    information - contains further information extracted from commonroad
%
% Example:
%    [statObs,dynObs,x0,goalSet,lanelets] = commonroad2cora('USA_US101-1_1_S-1', 'verbose', false);
%
%    figure
%    hold on
%    for i = 1:length(lanelets)
%       plot(lanelets{i},[1,2],'FaceColor',[.7 .7 .7],'Filled',true);
%    end
%    plot(goalSet{1}.set,[1,2],'r','EdgeColor','none','FaceAlpha',0.5,'Filled',true);
%    plot(x0.x,x0.y,'.g','MarkerSize',20);
%    for i = 1:length(dynObs)
%       plot(dynObs{i}.set,[1,2],'b','EdgeColor','none','Filled',true);
%    end
%
%    xlim([-32,35]);
%    ylim([-14,15]);
%    axis equal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: spaceex2cora

% Author:       Farah Atour, Niklas Kochdumper, Philipp Gassert
% Written:      28-April-2020
% Last update:  07-August-2020
% Last revision: ---

%------------- BEGIN CODE --------------%

tic;

% ----------------------------
% ----- Reading options ------
% ----------------------------

defaultVerbose = true;

p = inputParser;
addParameter(p, 'verbose', defaultVerbose, @(x) islogical(x) || (x == 1) || (x == 0) );

parse(p, varargin{:});

verbose = p.Results.verbose;

% temporary advisory

fprintf('-----------------------------\nConverting...\nPlease report errors to Niklas Kochdumper or Philipp Gassert under ga96man@mytum.de.\n');

% ----------------------------
% --- Reading the XML file ---
% ----------------------------
if ~contains(filename, '.xml')
    filename = strcat(filename,'.xml');
end
xDoc = xmlread(filename);

if xDoc.getElementsByTagName('commonRoad').getLength > 1
    throw(MException('Currently cannot handle more than one scenario!'));
end

% --- Retreive commonRoad version ---
commonroad_version = xDoc.getElementsByTagName('commonRoad').item(0).getAttribute('commonRoadVersion');

% -----------------------------------------------
% -- Extract all lanelets of the road network ---
% -----------------------------------------------

laneletList = xDoc.getElementsByTagName('lanelet');

% --- Preallocate for better memory and speed ---
lanelet_length = laneletList.getLength;
laneletFlag = lanelet_length > 0;
lanelets = struct('id', cell(lanelet_length,1), ...
    'leftBound', cell(lanelet_length,1), ...
    'rightBound', cell(lanelet_length,1));
pgon = polyshape();


for i = 0:(lanelet_length-1)
    % id
    if ~strcmp(laneletList.item(i).getAttribute('id'),'')
        lanelets(i+1).id = str2double(laneletList.item(i).getAttribute('id'));
    else
        continue
    end
    
    % left bound
    boundList = laneletList.item(i).getElementsByTagName('leftBound');
    leftbound_pointList = boundList.item(0).getElementsByTagName('point');
    lanelets(i+1).leftBound = zeros(2,leftbound_pointList.getLength());
    for j = 0:(leftbound_pointList.getLength()-1)
        lanelets(i+1).leftBound(1,j+1) = str2double(leftbound_pointList.item(j).getElementsByTagName('x').item(0).getTextContent);
        lanelets(i+1).leftBound(2,j+1) = str2double(leftbound_pointList.item(j).getElementsByTagName('y').item(0).getTextContent);
    end
    
    % right bound
    boundList = laneletList.item(i).getElementsByTagName('rightBound');
    rightbound_pointList = boundList.item(0).getElementsByTagName('point');
    lanelets(i+1).rightBound = zeros(2,rightbound_pointList.getLength());
    for j = 0:(rightbound_pointList.getLength()-1)
        lanelets(i+1).rightBound(1,j+1) = str2double(rightbound_pointList.item(j).getElementsByTagName('x').item(0).getTextContent);
        lanelets(i+1).rightBound(2,j+1) = str2double(rightbound_pointList.item(j).getElementsByTagName('y').item(0).getTextContent);
    end
    
    % vertices of the left bound lanelet
    x_left = lanelets(i+1).leftBound(1,:);
    y_left = lanelets(i+1).leftBound(2,:);
    
    % vertices of the right bound lanelet
    x_right = flip(lanelets(i+1).rightBound(1,:));
    y_right = flip(lanelets(i+1).rightBound(2,:));
    
    % creating a polygon of the lanelet
    pgon(i+1) = polyshape([x_left,x_right],[y_left,y_right],'Simplify',false);
    lanelet_poly(i+1).polygon = pgon(i+1).Vertices;
    
    % adding polygons to get one polygon of the network road
    if i == 0
        poly = pgon(i+1);
    else
        poly = union(poly,pgon(i+1));
    end
end

% -----------------------------------------------
% -- Extract all static and dynamic obstacles ---
% -----------------------------------------------

timeStep = str2double(xDoc.item(0).getAttribute('timeStepSize'));
dynamic_counter = 0;
static_counter  = 0;
dynamicList =  struct();
staticList =  struct();
trajectory_count = 0;
occupancy_count = 0;

% --- Extracting obstacles for 2020a version; dynamicObstacles/staticObstacles separate ---
if strcmp(commonroad_version, '2020a')
    dynamicObstaclesList = xDoc.getElementsByTagName('dynamicObstacle');
    dynamicObstacles_length = dynamicObstaclesList.getLength;
    staticObstaclesList = xDoc.getElementsByTagName('staticObstacle');
    staticObstacles_length = staticObstaclesList.getLength;
    
    % extracting dynamic obstacles
    if ~isempty(dynamicObstaclesList.item(0))
        dynamicObstacles(dynamicObstacles_length).id = [];
        for i = 0:(dynamicObstacles_length-1)
            % id
            dynamicObstacles(i+1).id = str2double(dynamicObstaclesList.item(i).getAttribute('id'));
            
            % initial state
            InitialStateList = dynamicObstaclesList.item(i).getElementsByTagName('initialState');
            ShapeList_initialstate = dynamicObstaclesList.item(i).getElementsByTagName('shape');
            positionList = InitialStateList.item(0).getElementsByTagName('position');
            orient = str2double(InitialStateList.item(0).getElementsByTagName('orientation').item(0).getTextContent);
            pointList_initialstate = positionList.item(0).getElementsByTagName('point');
            points_initialstate = points_fct(pointList_initialstate);
            x_initialstate = points_initialstate.x;
            y_initialstate = points_initialstate.y;
            % transforming any shape into a polygon
            polygon_initialstate = shape_fct(ShapeList_initialstate,x_initialstate,y_initialstate);
            temp = polyshape(polygon_initialstate,'Simplify',false);
            temp = rotate(temp,rad2deg(orient),[x_initialstate,y_initialstate]);
            initialstate_vertices = temp.Vertices;
            
            dynamic_counter = dynamic_counter + 1;
            % occupancySet
            occupancySetList = dynamicObstaclesList.item(i).getElementsByTagName('occupancySet');
            if ~isempty(occupancySetList.item(0))
                occupancy_count = occupancy_count + 1;
                occupancyList = occupancySetList.item(0).getElementsByTagName('occupancy');
                % inserting the rest
                x = x_initialstate;
                y = y_initialstate;
                dynamicObstacles(i+1).dynamic = occupancy_fct(occupancyList,x,y);
                dynamicList(dynamic_counter).dynamic = occupancy_fct(occupancyList,x,y);
                % inserting the initial state
                dynamicObstacles(i+1).dynamic(1).polygon = initialstate_vertices;
                dynamicList(dynamic_counter).dynamic(1).polygon = initialstate_vertices;
                dynamicObstacles(i+1).dynamic(1).time = 0;
                dynamicList(dynamic_counter).dynamic(1).time = 0;
            end
            
            % trajectory
            trajectoryList = dynamicObstaclesList.item(i).getElementsByTagName('trajectory');
            if ~isempty(trajectoryList.item(0))
                trajectory_count = trajectory_count + 1;
                trajectorystateList = trajectoryList.item(0).getElementsByTagName('state');
                % inserting the rest
                dynamicObstacles(i+1).dynamic = trajectory_fct(trajectorystateList,ShapeList_initialstate);
                dynamicList(dynamic_counter).dynamic = trajectory_fct(trajectorystateList,ShapeList_initialstate);
                % inserting the initial state
                dynamicObstacles(i+1).dynamic(1).polygon = initialstate_vertices;
                dynamicList(dynamic_counter).dynamic(1).polygon = initialstate_vertices;
                dynamicObstacles(i+1).dynamic(1).time = 0;
                dynamicList(dynamic_counter).dynamic(1).time = 0;
            end
        end
    end
    
    % extracting static obstacles
    if ~isempty(staticObstaclesList.item(0))
        staticObstacles(staticObstacles_length).id = [];
        for i = 0:(staticObstacles_length-1)
            % id
            staticObstacles(i+1).id = str2double(staticObstaclesList.item(i).getAttribute('id'));
            
            % initial state
            InitialStateList = staticObstaclesList.item(i).getElementsByTagName('initialState');
            ShapeList_initialstate = staticObstaclesList.item(i).getElementsByTagName('shape');
            positionList = InitialStateList.item(0).getElementsByTagName('position');
            orient = str2double(InitialStateList.item(0).getElementsByTagName('orientation').item(0).getTextContent);
            pointList_initialstate = positionList.item(0).getElementsByTagName('point');
            points_initialstate = points_fct(pointList_initialstate);
            x_initialstate = points_initialstate.x;
            y_initialstate = points_initialstate.y;
            % transforming any shape into a polygon
            polygon_initialstate = shape_fct(ShapeList_initialstate,x_initialstate,y_initialstate);
            temp = polyshape(polygon_initialstate,'Simplify',false);
            temp = rotate(temp,rad2deg(orient),[x_initialstate,y_initialstate]);
            
            static_counter = static_counter + 1;
            
            % obstacle is a polygon
            polygonList = ShapeList_initialstate.item(0).getElementsByTagName('polygon');
            
            for l = 0:(polygonList.getLength()-1)
                polygon_pointList = polygonList.item(l).getElementsByTagName('point');
                x = zeros(1,polygon_pointList.getLength());
                y = zeros(1,polygon_pointList.getLength());
                for k = 0:(polygon_pointList.getLength()-1)
                    x(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('x').item(0).getTextContent);
                    y(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('y').item(0).getTextContent);
                end
                temp = polyshape(x,y,'Simplify',false);
                if l == 0
                    shape_polygon = temp;
                else
                    shape_polygon = union(shape_polygon,temp);
                end
            end
            
            % obstacle is a rectangle
            rectangleList = ShapeList_initialstate.item(0).getElementsByTagName('rectangle');
            
            rectangle = polyshape();
            for l = 0:(rectangleList.getLength()-1)
                xCenter = str2double(rectangleList.item(l).getElementsByTagName('center').item(0).getElementsByTagName('x').item(0).getTextContent);
                xCenter = xCenter + str2double(InitialStateList.item(l).getElementsByTagName('position').item(0).getElementsByTagName('point').item(0).getElementsByTagName('x').item(0).getTextContent);
                yCenter = str2double(rectangleList.item(l).getElementsByTagName('center').item(0).getElementsByTagName('y').item(0).getTextContent);
                yCenter = yCenter + str2double(InitialStateList.item(l).getElementsByTagName('position').item(0).getElementsByTagName('point').item(0).getElementsByTagName('y').item(0).getTextContent);
                width = str2double(rectangleList.item(l).getElementsByTagName('length').item(0).getTextContent);
                height = str2double(rectangleList.item(l).getElementsByTagName('width').item(0).getTextContent);
                orient = str2double(rectangleList.item(l).getElementsByTagName('orientation').item(0).getTextContent);
                orient = orient + str2double(InitialStateList.item(l).getElementsByTagName('orientation').item(0).getElementsByTagName('exact').item(0).getTextContent);
                
                w = width/2;
                h = height/2;
                xLeftBottom = xCenter - cos(orient) * w + sin(orient) * h;
                yLeftBottom = yCenter - cos(orient) * h - sin(orient) * w;
                xRightBottom = xCenter + cos(orient) * w + sin(orient) * h;
                yRightBottom = yCenter - cos(orient) * h + sin(orient) * w;
                xLeftTop = xCenter - cos(orient) * w - sin(orient) * h;
                yLeftTop = yCenter + cos(orient) * h - sin(orient) * w;
                xRightTop = xCenter + cos(orient) * w - sin(orient) * h;
                yRightTop = yCenter + cos(orient) * h + sin(orient) * w;
                % convert rectangle to corner points list
                points = [xLeftBottom xLeftTop xRightTop xRightBottom; ...
                    yLeftBottom yLeftTop yRightTop yRightBottom]';
                % converting the rectangle into a polygon
                rectangle(l+1) = polyshape(points,'Simplify',false);
                % sum all the rectangles(polygons)
                if l == 0
                    shape_polygon = rectangle(l+1);
                else
                    shape_polygon = union(shape_polygon,rectangle(l+1));
                end
            end
            
            % obstacle is a circle
            circleList = ShapeList_initialstate.item(0).getElementsByTagName('circle');
            
            % "unit test" for dfinition of shape
            if (polygonList.getLength > 0) && (circleList.getLength > 0)
                disp('Error Warning: Obstacle is both circle and polygon! Cannot handle: Will only handle circle.')
            end
            
            circle = polyshape();
            for l = 0:(circleList.getLength()-1)
                xCenter = str2double(circleList.item(0).getElementsByTagName('center').item(0).getElementsByTagName('x').item(0).getTextContent);
                yCenter = str2double(circleList.item(0).getElementsByTagName('center').item(0).getElementsByTagName('y').item(0).getTextContent);
                radius = str2double(circleList.item(l).getElementsByTagName('radius').item(0).getTextContent);
                % convering circle into a polygon
                n = 50; %number of points
                theta = (0:n-1)*(2*pi/n);
                x = xCenter + radius*cos(theta);
                y = yCenter + radius*sin(theta);
                circle(l+1) = polyshape(x,y,'Simplify',false);
                % sum all the circles(polygons)
                if l == 0
                    shape_polygon = circle(l+1);
                else
                    shape_polygon = union(shape_polygon,circle(l+1));
                end
            end
            
            reg = regions(shape_polygon);
            
            for l = 1:size(reg,1)
                vert = reg(l,1).Vertices;
                staticObstacles(i+1).static(l).polygon = vert;
                staticList(static_counter).static(l).polygon = vert;
                staticObstacles(i+1).static(l).time = 0;
                staticList(static_counter).static(l).time = 0;
            end
        end
    end
    
    % --- Extracting obstacles from 2018b version ---
elseif strcmp(commonroad_version, '2018b')
    obstaclesList = xDoc.getElementsByTagName('obstacle');
    obstacles_length = obstaclesList.getLength;
    if ~isempty(obstaclesList.item(0))
        obstacles(obstacles_length).id = [];
        for i = 0:(obstacles_length-1)
            % id
            obstacles(i+1).id = str2double(obstaclesList.item(i).getAttribute('id'));
            
            % role
            roleList = obstaclesList.item(i).getElementsByTagName('role');
            role = char(roleList.item(0).getTextContent);
            
            % initial state
            InitialStateList = obstaclesList.item(i).getElementsByTagName('initialState');
            ShapeList_initialstate = obstaclesList.item(i).getElementsByTagName('shape');
            positionList = InitialStateList.item(0).getElementsByTagName('position');
            orient = str2double(InitialStateList.item(0).getElementsByTagName('orientation').item(0).getTextContent);
            pointList_initialstate = positionList.item(0).getElementsByTagName('point');
            points_initialstate = points_fct(pointList_initialstate);
            x_initialstate = points_initialstate.x;
            y_initialstate = points_initialstate.y;
            % transforming any shape into a polygon
            polygon_initialstate = shape_fct(ShapeList_initialstate,x_initialstate,y_initialstate);
            temp = polyshape(polygon_initialstate,'Simplify',false);
            temp = rotate(temp,rad2deg(orient),[x_initialstate,y_initialstate]);
            initialstate_vertices = temp.Vertices;
            
            if ~strcmp(role, 'static')
                dynamic_counter = dynamic_counter + 1;
                % occupancySet
                occupancySetList = obstaclesList.item(i).getElementsByTagName('occupancySet');
                if ~isempty(occupancySetList.item(0))
                    occupancy_count = occupancy_count + 1;
                    occupancyList = occupancySetList.item(0).getElementsByTagName('occupancy');
                    % inserting the rest
                    x = x_initialstate;
                    y = y_initialstate;
                    obstacles(i+1).dynamic = occupancy_fct(occupancyList,x,y);
                    dynamicList(dynamic_counter).dynamic = occupancy_fct(occupancyList,x,y);
                    % inserting the initial state
                    obstacles(i+1).dynamic(1).polygon = initialstate_vertices;
                    dynamicList(dynamic_counter).dynamic(1).polygon = initialstate_vertices;
                    obstacles(i+1).dynamic(1).time = 0;
                    dynamicList(dynamic_counter).dynamic(1).time = 0;
                end
                
                % trajectory
                trajectoryList = obstaclesList.item(i).getElementsByTagName('trajectory');
                if ~isempty(trajectoryList.item(0))
                    trajectory_count = trajectory_count + 1;                    
                    trajectorystateList = trajectoryList.item(0).getElementsByTagName('state');
                    % inserting the rest
                    obstacles(i+1).dynamic = trajectory_fct(trajectorystateList,ShapeList_initialstate);
                    dynamicList(dynamic_counter).dynamic = trajectory_fct(trajectorystateList,ShapeList_initialstate);
                    % inserting the initial state
                    obstacles(i+1).dynamic(1).polygon = initialstate_vertices;
                    dynamicList(dynamic_counter).dynamic(1).polygon = initialstate_vertices;
                    obstacles(i+1).dynamic(1).time = 0;
                    dynamicList(dynamic_counter).dynamic(1).time = 0;
                end
                
            else
                
                % static obstacle
                static_counter = static_counter + 1;
                
                % obstacle is a polygon
                polygonList = ShapeList_initialstate.item(0).getElementsByTagName('polygon');
                
                for l = 0:(polygonList.getLength()-1)
                    polygon_pointList = polygonList.item(l).getElementsByTagName('point');
                    x = zeros(1,polygon_pointList.getLength());
                    y = zeros(1,polygon_pointList.getLength());
                    for k = 0:(polygon_pointList.getLength()-1)
                        x(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('x').item(0).getTextContent);
                        y(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('y').item(0).getTextContent);
                    end
                    temp = polyshape(x,y,'Simplify',false);
                    if l == 0
                        shape_polygon = temp;
                    else
                        shape_polygon = union(shape_polygon,temp);
                    end
                end
                
                % obstacle is a rectangle
                rectangleList = ShapeList_initialstate.item(0).getElementsByTagName('rectangle');
                
                rectangle = polyshape();
                for l = 0:(rectangleList.getLength()-1)
                    xCenter = str2double(rectangleList.item(l).getElementsByTagName('center').item(0).getElementsByTagName('x').item(0).getTextContent);
                    xCenter = xCenter + str2double(InitialStateList.item(l).getElementsByTagName('position').item(0).getElementsByTagName('point').item(0).getElementsByTagName('x').item(0).getTextContent);
                    yCenter = str2double(rectangleList.item(l).getElementsByTagName('center').item(0).getElementsByTagName('y').item(0).getTextContent);
                    yCenter = yCenter + str2double(InitialStateList.item(l).getElementsByTagName('position').item(0).getElementsByTagName('point').item(0).getElementsByTagName('y').item(0).getTextContent);
                    width = str2double(rectangleList.item(l).getElementsByTagName('length').item(0).getTextContent);
                    height = str2double(rectangleList.item(l).getElementsByTagName('width').item(0).getTextContent);
                    orient = str2double(rectangleList.item(l).getElementsByTagName('orientation').item(0).getTextContent);
                    orient = orient + str2double(InitialStateList.item(l).getElementsByTagName('orientation').item(0).getElementsByTagName('exact').item(0).getTextContent);

                    w = width/2;
                    h = height/2;
                    xLeftBottom = xCenter - cos(orient) * w + sin(orient) * h;
                    yLeftBottom = yCenter - cos(orient) * h - sin(orient) * w;
                    xRightBottom = xCenter + cos(orient) * w + sin(orient) * h;
                    yRightBottom = yCenter - cos(orient) * h + sin(orient) * w;
                    xLeftTop = xCenter - cos(orient) * w - sin(orient) * h;
                    yLeftTop = yCenter + cos(orient) * h - sin(orient) * w;
                    xRightTop = xCenter + cos(orient) * w - sin(orient) * h;
                    yRightTop = yCenter + cos(orient) * h + sin(orient) * w;
                    % convert rectangle to corner points list
                    points = [xLeftBottom xLeftTop xRightTop xRightBottom; ...
                        yLeftBottom yLeftTop yRightTop yRightBottom]';
                    % converting the rectangle into a polygon
                    rectangle(l+1) = polyshape(points,'Simplify',false);
                    % sum all the rectangles(polygons)
                    if l == 0
                        shape_polygon = rectangle(l+1);
                    else
                        shape_polygon = union(shape_polygon,rectangle(l+1));
                    end
                end
                
                % obstacle is a circle
                circleList = ShapeList_initialstate.item(0).getElementsByTagName('circle');
                
                circle = polyshape();
                for l = 0:(circleList.getLength()-1)
                    xCenter = str2double(circleList.item(0).getElementsByTagName('center').item(0).getElementsByTagName('x').item(0).getTextContent);
                    yCenter = str2double(circleList.item(0).getElementsByTagName('center').item(0).getElementsByTagName('y').item(0).getTextContent);
                    radius = str2double(circleList.item(l).getElementsByTagName('radius').item(0).getTextContent);
                    % convering circle into a polygon
                    n = 50; %number of points
                    theta = (0:n-1)*(2*pi/n);
                    x = xCenter + radius*cos(theta);
                    y = yCenter + radius*sin(theta);
                    circle(l+1) = polyshape(x,y,'Simplify',false);
                    % sum all the circles(polygons)
                    if l == 0
                        shape_polygon = circle(l+1);
                    else
                        shape_polygon = union(shape_polygon,circle(l+1));
                    end
                end
                
                reg = regions(shape_polygon);
                
                for l = 1:size(reg,1)
                    vert = reg(l,1).Vertices;
                    obstacles(i+1).static(l).polygon = vert;
                    staticList(static_counter).static(l).polygon = vert;
                    obstacles(i+1).static(l).time = 0;
                    staticList(static_counter).static(l).time = 0;
                end
            end
        end
    end
    
else
    disp('Unidentified version (not 2018b or 2020a) when extracting obstacle information. No obstacles extracted.');
end

% -----------------------------------------------
% -- Extract all Ego Vehicles -------------------
% -----------------------------------------------

egoVehiclesList = xDoc.getElementsByTagName('planningProblem');
if ~isempty(egoVehiclesList.item(0))
    egoVehicles(egoVehiclesList.getLength()).initial = [];
    for i = 0:(egoVehiclesList.getLength()-1)
        
        % initial state
        InitialStateList = egoVehiclesList.item(i).getElementsByTagName('initialState');
        positionList = InitialStateList.item(0).getElementsByTagName('position');
        pointList_initialstate = positionList.item(0).getElementsByTagName('point');
        points_initialstate = points_fct(pointList_initialstate);
        egoVehicles(i+1).initial.x = points_initialstate.x;
        egoVehicles(i+1).initial.y = points_initialstate.y;
        
        timeList_initial = InitialStateList.item(0).getElementsByTagName('time');
        if ~isempty(timeList_initial.item(0).getElementsByTagName('exact').item(0))
            egoVehicles(i+1).initial.time = str2double(timeList_initial.item(0).getElementsByTagName('exact').item(0).getTextContent);
        elseif ~isempty(timeList_initial.item(0).getElementsByTagName('intervalStart').item(0))
            egoVehicles(i+1).initial.time(1) = str2double(timeList_initial.item(0).getElementsByTagName('intervalStart').item(0).getTextContent);
            egoVehicles(i+1).initial.time(2) = str2double(timeList_initial.item(0).getElementsByTagName('intervalEnd').item(0).getTextContent);
        end
        
        velocity_initial = InitialStateList.item(0).getElementsByTagName('velocity');
        egoVehicles.initial.velocity = str2double(velocity_initial.item(0).getElementsByTagName('exact').item(0).getTextContent);
        
        orientation_initial = InitialStateList.item(0).getElementsByTagName('orientation');
        egoVehicles.initial.orientation = str2double(orientation_initial.item(0).getElementsByTagName('exact').item(0).getTextContent);
        
        % goal region
        goalStateList = egoVehiclesList.item(i).getElementsByTagName('goalState');
        
        % position
        positionList = goalStateList.item(0).getElementsByTagName('position');
        if ~isempty(positionList.item(0))
            % point
            pointList_goal = positionList.item(0).getElementsByTagName('point');
            if ~isempty(pointList_goal.item(0))
                points_goal = points_fct(pointList_goal);
                egoVehicles(i+1).goalRegion.x = points_goal.x;
                egoVehicles(i+1).goalRegion.y = points_goal.y;
            end
            % rectangle
            rectangleList_goal = positionList.item(0).getElementsByTagName('rectangle');
            if ~isempty(rectangleList_goal.item(0))
                wid = str2double(rectangleList_goal.item(0).getElementsByTagName('length').item(0).getTextContent);
                len = str2double(rectangleList_goal.item(0).getElementsByTagName('width').item(0).getTextContent);
                orient = str2double(rectangleList_goal.item(0).getElementsByTagName('orientation').item(0).getTextContent);
                centerpoints =rectangleList_goal.item(0).getElementsByTagName('center');
                x = str2double(centerpoints.item(0).getElementsByTagName('x').item(0).getTextContent);
                y = str2double(centerpoints.item(0).getElementsByTagName('y').item(0).getTextContent);
                x_left = x - wid / 2;
                y_bottom = y - len / 2;
                % convert rectangle to corner points list
                points_rec = [x_left x_left x_left + wid x_left + wid; ...
                    y_bottom y_bottom + len y_bottom + len y_bottom]';
                % converting the rectangle into a polygon
                temp = polyshape(points_rec,'Simplify',false);
                temp = rotate(temp,rad2deg(orient),[x,y]);
                egoVehicles(i+1).goalRegion.polygon.Poly = temp.Vertices;
            end
            % lanelet
            laneletList = positionList.item(0).getElementsByTagName('lanelet');
            lanelet_length = laneletList.getLength;
            for m = 0:(lanelet_length-1)
                ref = str2double(laneletList.item(m).getAttribute('ref'));
                list = [lanelets.id];
                idx = list == ref;
                egoVehicles(i+1).goalRegion.polygon(m+1).Poly=  pgon(idx).Vertices;
            end
        end
        
        % time
        timeList_goal = goalStateList.item(0).getElementsByTagName('time');
        if ~isempty(timeList_goal.item(0).getElementsByTagName('exact').item(0))
            egoVehicles(i+1).goalRegion.time = str2double(timeList_goal.item(0).getElementsByTagName('exact').item(0).getTextContent);
        elseif ~isempty(timeList_goal.item(0).getElementsByTagName('intervalStart').item(0))
            egoVehicles(i+1).goalRegion.time(1) = str2double(timeList_goal.item(0).getElementsByTagName('intervalStart').item(0).getTextContent);
            egoVehicles(i+1).goalRegion.time(2) = str2double(timeList_goal.item(0).getElementsByTagName('intervalEnd').item(0).getTextContent);
        end
        
        % velocity
        velList_goal = goalStateList.item(0).getElementsByTagName('velocity');
        if ~isempty(velList_goal.item(0))
            if ~isempty(velList_goal.item(0).getElementsByTagName('exact').item(0))
                egoVehicles(i+1).goalRegion.velocity = str2double(velList_goal.item(0).getElementsByTagName('exact').item(0).getTextContent);
            elseif ~isempty(velList_goal.item(0).getElementsByTagName('intervalStart').item(0))
                egoVehicles(i+1).goalRegion.velocity(1) = str2double(velList_goal.item(0).getElementsByTagName('intervalStart').item(0).getTextContent);
                egoVehicles(i+1).goalRegion.velocity(2) = str2double(velList_goal.item(0).getElementsByTagName('intervalEnd').item(0).getTextContent);
            end
        end
        
        % orientation
        orientList_goal = goalStateList.item(0).getElementsByTagName('orientation');
        if ~isempty(orientList_goal.item(0))
            if ~isempty(orientList_goal.item(0).getElementsByTagName('exact').item(0))
                egoVehicles(i+1).goalRegion.orientation = str2double(orientList_goal.item(0).getElementsByTagName('exact').item(0).getTextContent);
            elseif ~isempty(orientList_goal.item(0).getElementsByTagName('intervalStart').item(0))
                egoVehicles(i+1).goalRegion.orientation(1) = str2double(orientList_goal.item(0).getElementsByTagName('intervalStart').item(0).getTextContent);
                egoVehicles(i+1).goalRegion.orientation(2) = str2double(orientList_goal.item(0).getElementsByTagName('intervalEnd').item(0).getTextContent);
            end
        end
    end
else
    egoVehicles = [];
end

% -----------------------------------------------
% -- Extract further information ----------------
% -----------------------------------------------

if strcmp(commonroad_version, '2020a')
    information = struct();
    
    % tags
    tags_all = xDoc.getElementsByTagName('scenarioTags').item(0).getElementsByTagName('*');
    tags = cell(length(tags_all));
    
    for i = 0:length(tags_all)
        tags{i+1} = string(tags_all.item(i).getTagName);
    end
    information.tags = tags;
    
    % location
    information.location.geoNameId = str2double(xDoc.getElementsByTagName('location').item(0).getElementsByTagName('geoNameId').item(0).getTextContent);
    information.location.gpsLatitude = str2double(xDoc.getElementsByTagName('location').item(0).getElementsByTagName('gpsLatitude').item(0).getTextContent);
    information.location.gpsLongitude = str2double(xDoc.getElementsByTagName('location').item(0).getElementsByTagName('gpsLongitude').item(0).getTextContent);
    
end

% -----------------------------------------------
% -- Transform to required output format --------
% -----------------------------------------------

w = warning();
warning('off');

% Static Obstacles
statObs = [];

% --- Retreive information if obstacles exist ---
if static_counter
    for i = 1:length(staticList)
        for j = 1:length(staticList(i).static)
            % get set
            points = staticList(i).static(j).polygon;
            statObs{end+1,1} = polygon(points(:,1),points(:,2));
        end
    end
end

% Dynamic Obstacles
dynObs = [];

% --- Retreive information if obstacles exist ---
if dynamic_counter
    for i = 1:length(dynamicList)
        for j = 1:length(dynamicList(i).dynamic)
            % get set
            points = dynamicList(i).dynamic(j).polygon;
            dynObs{end+1,1}.set = polygon(points(:,1),points(:,2));
            % get time interval;
            t = dynamicList(i).dynamic(j).time*timeStep;
            if length(t) == 1
                dynObs{end,1}.time = interval(t);
            else
                dynObs{end,1}.time = interval(t(1),t(2));
            end
        end
    end
end

% ---Retreive planning problem if it exists ---
if ~isempty(egoVehicles)
    
    % Initial Point
    x0 = egoVehicles.initial;
    x0.time = x0.time*timeStep;
    
    
    % Goal Region
    t = egoVehicles.goalRegion.time*timeStep;
    if length(t) == 1
        time = interval(t);
    else
        time = interval(t(1),t(2));
    end
    
    if isfield(egoVehicles.goalRegion,'velocity')
        v = egoVehicles.goalRegion.velocity;
        if length(v) == 1
            velocity = interval(v);
        else
            velocity = interval(v(1),v(2));
        end
    else
        velocity = [];
    end
    
    if isfield(egoVehicles.goalRegion,'orientation')
        o = egoVehicles.goalRegion.orientation;
        if length(o) == 1
            orientation = interval(o);
        else
            orientation = interval(o(1),o(2));
        end
    else
        orientation = [];
    end
    
    if isfield(egoVehicles.goalRegion,'polygon')
        goalSet = cell(length(egoVehicles.goalRegion.polygon),1);
        
        for i = 1:length(egoVehicles.goalRegion.polygon)
            
            % get goal set set
            points = egoVehicles.goalRegion.polygon(i).Poly;
            goalSet{i}.set = polygon(points(:,1),points(:,2));
            
            % store time, orientation, and velocity
            goalSet{i}.time = time;
            goalSet{i}.orientation = orientation;
            goalSet{i}.velocity = velocity;
        end
    else
        goalSet{1}.time = egoVehicles.goalRegion.time;
        goalSet{1}.set = [];
    end
else
    % Return empty arrays if no planning problem exists
    disp('No planning problem found. Will return empty vectors x0, goalSet');
    x0 = [];
    goalSet = [];
end

% Lanelets (if existent)

if laneletFlag
    lanelets = cell(length(lanelet_poly),1);
    
    for i = 1:length(lanelet_poly)
        points = lanelet_poly(i).polygon;
        lanelets{i} = polygon(points(:,1),points(:,2));
    end
else
    lanelets = {};
end

warning(w);

converterTime = toc;

if (occupancy_count + trajectory_count ~= dynamic_counter)
    fprintf('Warning: During conversion a mismatch for the number of obstacles was encountered, i.e. number of obstacles by trajectory + number of obstacles by occupancy set unequal total number of obstacles. Your scenario may contain dynamic obstacles without dynamics.\nIf not, this is a converter problem; contact Niklas Kochdumper or Philipp Gassert under ga96man@mytum.de.\n');
end

if verbose
    fprintf(strcat('-----------------------------\nSuccessfully converted\n', filename,...
        ' with\n%d lanelets,\n%d static obstacles\n%d dynamic obstacles by trajectory\n%d dynamic obstacles by occupancy set\n%d planning problems in\n%d seconds.\n-----------------------------\n'),...
        length(lanelets), static_counter, trajectory_count, occupancy_count, egoVehiclesList.getLength, converterTime);
end

end


% Auxiliary Functions -----------------------------------------------------

function vertices = shape_fct(shapeList,x,y, orientList)
% central position of all shapes
xCenter = x;
yCenter = y;
%rectangle
rectangleList = shapeList.item(0).getElementsByTagName('rectangle');
%for preallocation
rectangle = polyshape();
for l = 0:(rectangleList.getLength()-1)
    width = str2double(rectangleList.item(l).getElementsByTagName('length').item(0).getTextContent);
    height = str2double(rectangleList.item(l).getElementsByTagName('width').item(0).getTextContent);
    if exist('orientList', 'var')
        orient = str2double(orientList.item(l).getTextContent);
    else
        orient = 0;
    end
%     % finding specific points
%     xLeft = xCenter - cos(width / 2);
%     yBottom = yCenter - sin(height / 2);
%     % convert rectangle to corner points list
%     points = [xLeft xLeft xLeft+width xLeft+width; ...
%         yBottom yBottom+height yBottom+height yBottom]';
    
    w = width/2;
    h = height/2;
    xLeftBottom = xCenter - cos(orient) * w + sin(orient) * h;
    yLeftBottom = yCenter - cos(orient) * h - sin(orient) * w;
    xRightBottom = xCenter + cos(orient) * w + sin(orient) * h;
    yRightBottom = yCenter - cos(orient) * h + sin(orient) * w;
    xLeftTop = xCenter - cos(orient) * w - sin(orient) * h;
    yLeftTop = yCenter + cos(orient) * h - sin(orient) * w;
    xRightTop = xCenter + cos(orient) * w - sin(orient) * h;
    yRightTop = yCenter + cos(orient) * h + sin(orient) * w;
    % convert rectangle to corner points list
    points = [xLeftBottom xLeftTop xRightTop xRightBottom; ...
        yLeftBottom yLeftTop yRightTop yRightBottom]';
    % converting the rectangle into a polygon
    rectangle(l+1) = polyshape(points,'Simplify',false);
    % sum all the rectangles(polygons)
    if l == 0
        shape_polygon = rectangle(l+1);
    else
        shape_polygon = union(shape_polygon,rectangle(l+1));
    end
end
%circle
circleList = shapeList.item(0).getElementsByTagName('circle');
%for preallocation
circle = polyshape();
for l = 0:(circleList.getLength()-1)
    radius = str2double(circleList.item(l).getElementsByTagName('radius').item(0).getTextContent);
    % convering circle into a polygon
    n = 50; %number of points
    theta = (0:n-1)*(2*pi/n);
    x = xCenter + radius*cos(theta);
    y = yCenter + radius*sin(theta);
    circle(l+1) = polyshape(x,y,'Simplify',false);
    % sum all the circles(polygons)
    if l == 0
        shape_polygon = circle(l+1);
    else
        shape_polygon = union(shape_polygon,circle(l+1));
    end
end
%polygon
polygonList = shapeList.item(0).getElementsByTagName('polygon');
%for preallocation
pol = polyshape();
for l = 0:(polygonList.getLength()-1)
    polygon_pointList = polygonList.item(l).getElementsByTagName('point');
    x = zeros(1,polygon_pointList.getLength());
    y = zeros(1,polygon_pointList.getLength());
    for k = 0:(polygon_pointList.getLength()-1)
        x(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('x').item(0).getTextContent);
        y(k+1) = str2double(polygon_pointList.item(k).getElementsByTagName('y').item(0).getTextContent);
    end
    pol(l+1) = polyshape(x,y,'Simplify',false);
    if l == 0
        shape_polygon = pol(l+1);
    else
        shape_polygon = union(shape_polygon,pol(l+1));
    end
    
end
vertices =  shape_polygon.Vertices;
end

function occupancy = occupancy_fct(occupancyList,x,y)
occupancy(occupancyList.getLength).polygon = [];
for l = 0:(occupancyList.getLength()-1)
    % shape
    shapeList_occupancy = occupancyList.item(l).getElementsByTagName('shape');
    if ~isempty(shapeList_occupancy.item(0))
        occupancy(l+2).polygon = shape_fct(shapeList_occupancy,x,y);
    end
    % time
    timeList = occupancyList.item(l).getElementsByTagName('time');
    if ~isempty(timeList.item(0).getElementsByTagName('exact').item(0))
        occupancy(l+2).time = str2double(timeList.item(0).getElementsByTagName('exact').item(0).getTextContent);
    elseif ~isempty(timeList.item(0).getElementsByTagName('intervalStart').item(0))
        occupancy(l+2).time(1) = str2double(timeList.item(0).getElementsByTagName('intervalStart').item(0).getTextContent);
        occupancy(l+2).time(2) = str2double(timeList.item(0).getElementsByTagName('intervalEnd').item(0).getTextContent);
    end
end

end

function trajectory = trajectory_fct(trajectorystateList,ShapeList_initialstate)
trajectory(trajectorystateList.getLength).polygon = [];
for l = 0:(trajectorystateList.getLength()-1)
    % shape
    positionList_trajectory = trajectorystateList.item(l).getElementsByTagName('position');
    pointList = positionList_trajectory.item(0).getElementsByTagName('point');
    position = points_fct(pointList);
    orientList = trajectorystateList.item(l).getElementsByTagName('orientation');
    trajectory(l+2).polygon = shape_fct(ShapeList_initialstate,position.x,position.y, orientList);
    
    % time
    timeList = trajectorystateList.item(l).getElementsByTagName('time');
    if ~isempty(timeList.item(0).getElementsByTagName('exact').item(0))
        trajectory(l+2).time = str2double(timeList.item(0).getElementsByTagName('exact').item(0).getTextContent);
    end
end
end

function points = points_fct(pointList)

if ~isempty(pointList.item(0).getElementsByTagName('x').item(0))
    points.x = str2double(pointList.item(0).getElementsByTagName('x').item(0).getTextContent);
    points.y = str2double(pointList.item(0).getElementsByTagName('y').item(0).getTextContent);
else
    points.x = 0;
    points.y = 0;
end
end

