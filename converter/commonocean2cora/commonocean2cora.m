function [statObs,dynObs,x0,goalSet,waterways,shallows,information] = commonocean2cora(filename,varargin)
% commonocean2cora - convert CommonOcean XML-file scenario to CORA model
%
% Syntax:
%    [statObs,dynObs,x0,goalSet,waterways,shallows,information] = ...
%       commonocean2cora(filename);
%    [statObs,dynObs,x0,goalSet,waterways,shallows,information] = ...
%       commonocean2cora(filename,verbose);
%
% Inputs:
%    filename - name of the CommonOcean XML-file scenario (specified as string)
%    verbose - (optional) true/false whether conversion information should
%                         be printed to the terminal. Default is true.
%
% Outputs:
%    statObs - cell-array storing the static obstacles
%    dynObs - cell-array storing the set and the time interval for the
%             dynamic obstacles
%    x0 - initial point for the planning problem
%    goalSet - cell-array storing the goal sets and the corresponding time
%    waterways - cell-array storing the sets for all waterways
%    shallows - cell-array storing the shallows
%    information - contains further information extracted from commonocean
%
% Example:
%    [statObs,dynObs,x0,goalSet,waterways,shallows] = ...
%        commonocean2cora('ZAM_Tutorial-2_1_T-1',false);
% 
%    figure; hold on;
%    for i = 1:length(waterways)
%        plot(waterways{i}.set,'FaceColor',[.7 .7 .7]);
%    end
%    for i = 1:length(shallows)
%        plot(shallows{i}.set,'FaceColor','cyan');
%    end 
%    plot(goalSet{1}.set,[1,2],'FaceColor','r','FaceAlpha',0.5);
%    plot(x0.x,x0.y,'.g','MarkerSize',20);
%    for i = 1:length(dynObs)
%        plot(dynObs{i}.set,[1,2],'FaceColor','b','FaceAlpha',0.5);
%    end
%    for i = 1:length(statObs)
%        plot(statObs{i}.set,'FaceColor','y');
%    end
% 
%    xlim([-600,600]); ylim([-600,600]); axis equal;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Authors:       Farah Atour, Niklas Kochdumper, Philipp Gassert, Leni Rohe, Hanna Krasowski, Bruno Maione
% Written:       28-April-2020
% Last update:   16-February-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

conTimerStart = tic;

% set default values
verbose = setDefaultValues({true},varargin);

% check input arguments
inputArgsCheck({{filename,'att',{'char','string'}},...
                {verbose,'att','logical'}});

% ----------------------------
% --- Reading the XML file ---
% ----------------------------

% add '.xml' to the filename if missing
if ~contains(filename, '.xml')
    filename = strcat(filename,'.xml');
end
xDoc = xmlread(filename);

% assert number of scenarios
assert(xDoc.getElementsByTagName('commonOcean').getLength == 1, 'Cannot handle more than one scenario!');

% --- Retreive commonOcean version ---
commonocean_version = xDoc.getElementsByTagName('commonOcean').item(0).getAttribute('commonOceanVersion');

% -----------------------------------------------
% -- Extract all static and dynamic obstacles ---
% -----------------------------------------------

timeStep = str2double(xDoc.item(0).getAttribute('timeStepSize'));
dynamic_counter = 0;         % counts number of dynamic obstacles
static_counter  = 0;         % counts number of static obstacles
shallows_counter = 0;        % counts number of shallows
total_dynamic_counter = 0;   % counts total number of all (time-independent) dynObs
statObs = [];
dynObs = [];
shallows = [];

w = warning();
warning('off');

% --- Extracting obstacles for 2022a version ---
if strcmp(commonocean_version, '2022a')
    dynamicObstaclesList = xDoc.getElementsByTagName('dynamicObstacle');
    dynamicObstacles_length = dynamicObstaclesList.getLength;
    staticObstaclesList = xDoc.getElementsByTagName('staticObstacle');
    staticObstacles_length = staticObstaclesList.getLength;
    dynObsMat = [];
   
    counterBadOnes = 0;
    for i = 0:(dynamicObstacles_length-1)
        dynamic_counter = dynamic_counter + 1;
        
        % extracting dynamic obstacles
        dynamicBuffer = aux_dynamic_fct(dynamicObstaclesList, i);
        total_dynamic_counter = total_dynamic_counter + length(dynamicBuffer);
        
        for j = 1:length(dynamicBuffer)
            % get set
            points = dynamicBuffer(j).polygon;
            dynObs{end+1,1}.set = polygon(points(:,1),points(:,2));
            
            % get time interval;
            t = dynamicBuffer(j).time*timeStep;
            if length(t) == 1
                dynObs{end,1}.time = interval(t);
            else
                dynObs{end,1}.time = interval(t(1),t(2));
            end
        end

        % extracting dynamic obstracles matrix
        obsmat = aux_dynamic_mat(dynamicObstaclesList , timeStep, i);
        
        dynObsMat = [dynObsMat; obsmat];            
    end 

    % extracting static obstacles
    for i = 0:(staticObstacles_length-1)
        
        % initial state
        ShapeList_initialstate = aux_getXList(staticObstaclesList, {'shape'}, i);
        stateList_initialstate = aux_getXList(staticObstaclesList, {'initialState'}, i);
        
        % transforming any shape into a polyshape
        initialstate_vertices = aux_shape_fct(ShapeList_initialstate, stateList_initialstate);
        initialstate_polyshape = polyshape(initialstate_vertices); 
        
        % split into regions to obtain single-regions obstacles
        reg = regions(initialstate_polyshape);
        
        if length(reg) > 1
            fprintf('More than one region in static obstacle.\n')
        end
        
        for l = 1:size(reg,1)
            static_counter = static_counter + 1;
            vert = reg(l,1).Vertices;
            statObs{end+1,1} = polygon(vert(:,1),vert(:,2));
        end
    end
else
    fprintf('Unidentified version (not 2022a) when extracting obstacle information. No obstacles extracted.');
end

% -----------------------------------------------
% -- Extract all waterways of the Waters Network ---
% -----------------------------------------------

watersList = xDoc.getElementsByTagName('waterway');

% --- Preallocate for better memory and speed ---
waters_length = watersList.getLength();
generalWatersList = struct(...
    'id', cell(waters_length,1), ...
    'leftBound', cell(waters_length,1), ...
    'rightBound', cell(waters_length,1));
waterways = cell(waters_length,logical(waters_length));

for i = 0:(waters_length-1)
    
    % retrieve id and set to -1 for references from goal region specification
    if ~strcmp(watersList.item(i).getAttribute('id'),'')
        generalWatersList(i+1).id = str2double(watersList.item(i).getAttribute('id'));
    else
        generalWatersList(i+1).id = -1;
        continue
    end
    
    % left bound
    leftbound_pointList = aux_getXList(watersList, {'leftBound', 'point'}, i);
    generalWatersList(i+1).leftBound = zeros(2,leftbound_pointList.getLength());
    
    for j = 0:(leftbound_pointList.getLength()-1)
        generalWatersList(i+1).leftBound(1,j+1) = aux_getXEl(leftbound_pointList, {'x'}, j);
        generalWatersList(i+1).leftBound(2,j+1) = aux_getXEl(leftbound_pointList, {'y'}, j);
    end
    
    % right bound
    rightbound_pointList = aux_getXList(watersList, {'rightBound', 'point'}, i);
    generalWatersList(i+1).rightBound = zeros(2,rightbound_pointList.getLength());
    
    for j = 0:(rightbound_pointList.getLength()-1)
        generalWatersList(i+1).rightBound(1,j+1) = aux_getXEl(rightbound_pointList, {'x'}, j);
        generalWatersList(i+1).rightBound(2,j+1) = aux_getXEl(rightbound_pointList, {'y'}, j);
    end
    
    % vertices of the left bound lanelet (eases handling below)
    x_left = generalWatersList(i+1).leftBound(1,:);
    y_left = generalWatersList(i+1).leftBound(2,:);
    
    % vertices of the right bound lanelet (eases handling below)
    x_right = generalWatersList(i+1).rightBound(1,:);
    y_right = generalWatersList(i+1).rightBound(2,:);
    
    % insertion into output cell array
    waterways{i+1} = polygon([x_left, flip(x_right)], [y_left, flip(y_right)]);
    if waterways{i+1}.set.NumRegions > 1
        fprintf(strcat('Warning: Waterways consists of %d regions!',...
            ' Assuming imprecision contained in pointlist input file.',...
            ' Continuing without alteration.\n'), waterways{i+1}.set.NumRegions);
    end
end
   
% delete non-waters-constituting waters (i.e. waters references in goal
% regions that were picked up by getElementsByTagName)
for i = waters_length:-1:1
    if generalWatersList(i).id == -1
        generalWatersList(i) = [];
        waterways(i) = [];
    end
end

% -----------------------------------------------
% -- Extract all shallows of the Waters Network ---
% -----------------------------------------------

shallowsList = xDoc.getElementsByTagName('shallow');
shallows_length = shallowsList.getLength();

for i = 0:(shallows_length-1)

    % initial state
    ShapeList_initialstate = aux_getXList(shallowsList, {'shape'}, i);
    
    % transforming any shape into a polyshape
    initialstate_vertices = aux_shape_fct(ShapeList_initialstate);
    initialstate_polyshape = polyshape(initialstate_vertices); 

    % split into regions to obtain single-regions obstacles
    reg = regions(initialstate_polyshape);

    if length(reg) > 1
        fprintf('Warning: More than one region in shallow.\n')
    end

    for l = 1:size(reg,1)
        shallows_counter = shallows_counter + 1;
        vert = reg(l,1).Vertices;
        shallows{end+1,1} = polygon(vert(:,1),vert(:,2));
    end
end

    
% -----------------------------------------------
% -- Extract all Ego Vehicles -------------------
% -----------------------------------------------

num_goalRegions = 0;
num_planningProblems = 0;

egoVehiclesList = xDoc.getElementsByTagName('planningProblem');
if ~isempty(egoVehiclesList.item(0))
    num_planningProblems = egoVehiclesList.getLength();
    x0(num_planningProblems).x = [];
    goalSet = cell(1,num_planningProblems);
    
    for i = 0:(num_planningProblems-1)
        num_goalRegions(i+1) = 0;
        
        % initial state
        InitialStateList = egoVehiclesList.item(i).getElementsByTagName('initialState');
        x0(i+1).x = aux_getXEl(egoVehiclesList, {'initialState', 'position', 'point', 'x'}, i);
        x0(i+1).y = aux_getXEl(egoVehiclesList, {'initialState', 'position', 'point', 'y'}, i);
        
        % retrieval of x0-struct information
        timeList_initial = aux_getXList(InitialStateList, {'time'});
        x0(i+1).time = aux_interval_or_exact(timeList_initial);
        x0(i+1).velocity = aux_getXEl(InitialStateList, {'velocity', 'exact'});
        x0(i+1).orientation = aux_getXEl(InitialStateList, {'orientation', 'exact'});
        
        % goal region
        goalStateList = aux_getXList(egoVehiclesList, {'goalState'}, i);
        for k = 0:(goalStateList.getLength()-1)
            num_goalRegions(i+1) = num_goalRegions(i+1) + 1;
            
            % get position/goal space
            positionList = aux_getXList(goalStateList, {'position'}, k);
            if positionList.getLength() >= 2
                throw(CORAerror('CORA:converterIssue',...
                    'Violation of standard: More than one position for goalState.'));
            end
            if ~isempty(positionList.item(0))
                
                % waterways
                watersList = positionList.item(0).getElementsByTagName('waterway');
                waters_length = watersList.getLength;
                polyshapeBuffer = polyshape();
                for m = 0:(waters_length-1)
                    ref = str2double(watersList.item(m).getAttribute('ref'));
                    idx = [generalWatersList.id] == ref;
                    polyshapeBuffer = union(polyshapeBuffer,...
                        polyshape([generalWatersList(idx).leftBound, flip(generalWatersList(idx).rightBound, 2)]'));
                end
                
                % region by lanelets if found above
                if waters_length
                    goalSet{num_goalRegions(i+1),i+1}.set =...
                        polygon(polyshapeBuffer.Vertices(:,1), polyshapeBuffer.Vertices(:,2));
                
                % region by shape
                else
                    verticesBuffer = aux_shape_fct(positionList);
                    goalSet{num_goalRegions(i+1),i+1}.set =...
                        polygon(verticesBuffer(:,1), verticesBuffer(:,2));
                end
            else
                % insert empty set if no goal set found
                goalSet{num_goalRegions(i+1),i+1}.set = [];
            end
            
            % time
            timeList_goal = aux_getXList(goalStateList, {'time'}, k);
            goalSet{i+1}.time = aux_interval_or_exact(timeList_goal)*timeStep;
            
            % velocity
            velList_goal = aux_getXList(goalStateList, {'velocity'}, k);
            if ~isempty(velList_goal.item(0))
                goalSet{i+1}.velocity = aux_interval_or_exact(velList_goal);
            else
                goalSet{i+1}.velocity = [];
            end
            
            % orientation (Note: Additional measures have to be taken to
            % avoid confusion between shape orientation and state
            % orientation in xml
            orientList_goal = aux_getXList(goalStateList, {'orientation'}, k);
            if ~isempty(orientList_goal.item(0))
                for l = 0:(orientList_goal.getLength()-1)
                    if strcmp('goalState', orientList_goal.item(l).getParentNode.getNodeName)
                        goalSet{i+1}.orientation = aux_interval_or_exact(orientList_goal, l);
                        break
                    end
                end
            else
                goalSet{i+1}.orientation = [];
            end
        end
    end
else
    % Return empty arrays if no planning problem exists
    if verbose
        %fprintf('No planning problem found. Will return empty vectors x0, goalSet');
    end
    x0 = [];
    goalSet = [];
end

% -----------------------------------------------
% -- Extract further information ----------------
% -----------------------------------------------

if strcmp(commonocean_version, '2022a')
    information = struct();
    information.BadRes = counterBadOnes;
    information.timeStep = timeStep;
    try
        
        % tags
        tags_all = xDoc.getElementsByTagName('scenarioTags').item(0).getElementsByTagName('*');
        tags = cell(length(tags_all));
        
        for i = 0:(length(tags_all)-1)
            tags{i+1} = string(tags_all.item(i).getTagName);
        end
        information.tags = tags;

        % location
        information.location.geoNameId = str2double(xDoc.getElementsByTagName('location').item(0)...
            .getElementsByTagName('geoNameId').item(0).getTextContent);
        information.location.gpsLatitude = str2double(xDoc.getElementsByTagName('location').item(0)...
            .getElementsByTagName('gpsLatitude').item(0).getTextContent);
        information.location.gpsLongitude = str2double(xDoc.getElementsByTagName('location').item(0)...
            .getElementsByTagName('gpsLongitude').item(0).getTextContent);

    end
end

% -----------------------------------------------
% -- Transform to required output format --------
% -----------------------------------------------

% sanity checks
if length(statObs) ~= static_counter
    throw(CORAerror('CORA:converterIssue',...
        'Counter mismatch for static obstacles.'));
elseif length(dynObs) ~= total_dynamic_counter
    throw(CORAerror('CORA:converterIssue',...
        'Counter mismatch for dynamic obstacles.'));
end

% --- Check for multiple Planning Problems and reduce to one ---
if num_planningProblems > 1
    if verbose
        fprintf(strcat('There are several planning problems found. However, only one',...
            'can be handled by the converter, thus only the first one will be regarded.\n'));
    end
    x0 = x0(1);
    goalSet = goalSet(:,1);
end

warning(w);
converterTime = toc(conTimerStart);

if verbose
    fprintf(strcat('--------------\nSuccessfully converted\n', filename,...
        ' with\n%d waterways,\n%d shallows,\n',...
        '%d static obstacles,\n%d dynamic obstacles,\n',...
        '%d planning problem with\n%d goal regions in\n%.2f seconds.\n--------------\n'),...
        length(waterways),shallows_counter, static_counter, dynamic_counter, egoVehiclesList.getLength, num_goalRegions, converterTime);
end

end


% Auxiliary functions -----------------------------------------------------

function vertices = aux_shape_fct(shapeList, varargin)
% This function returns vertices for the combined polyshape which is the
% union of all shapes contained in shapeList, possibly combined with
% position and orientation from the optional parameter stateXList. You can
% choose which states from stateXList with index idx are relevant by
% passing the optional parameter idx.

[stateXList,idx] = setDefaultValues({[],0},varargin);

if isempty(stateXList)
    orientation_state = 0;
    x_state = 0;
    y_state = 0;
else
    if stateXList.getLength() < idx + 1
        throw(CORAerror('CORA:converterIssue',...
            sprintf('Expected state list of min length %i; got length %i.\n', idx, stateXList.getLength())));
    end
    orientation_state = aux_getXEl(stateXList, {'orientation', 'exact'}, idx);
    x_state = aux_getXEl(stateXList, {'position', 'point', 'x'}, idx);
    y_state = aux_getXEl(stateXList, {'position', 'point', 'y'}, idx);
end

% initialize empty shape
shape_polygon = polyshape();

% rectangle
rectangleList = aux_getXList(shapeList, {'rectangle'});
for l = 0:(rectangleList.getLength()-1)
    len = aux_getXEl(rectangleList, {'length'}, l);
    wid = aux_getXEl(rectangleList, {'width'}, l);
    
    if ~isempty(rectangleList.item(l).getElementsByTagName('center').item(0))
        x_center = aux_getXEl(rectangleList, {'center', 'x'}) + x_state;
        y_center = aux_getXEl(rectangleList, {'center', 'y'}) + y_state;
    else
        x_center = x_state;
        y_center = y_state;
    end
    
    orientationXList = aux_getXList(rectangleList, {'orientation'}, l);
    if ~isempty(orientationXList.item(0))
        if isempty(orientationXList.item(0).getElementsByTagName('exact').item(0))
            orientation = aux_getXEl(orientationXList) + orientation_state;
        else
            orientation = aux_getXEl(orientationXList, {'exact'}) + orientation_state;
        end
    else
        orientation = orientation_state;
    end
    
    points = [-len, -wid; len, -wid; len, wid; -len, wid]/2;
    points = points * [cos(orientation), sin(orientation); -sin(orientation), cos(orientation)] + [x_center, y_center];
    
    % add all the rectangles(polygons)
    shape_polygon = union(shape_polygon, polyshape(points,'Simplify',false));
end

% circle
circleList = aux_getXList(shapeList, {'circle'});
for l = 0:(circleList.getLength()-1)
    radius = aux_getXEl(circleList, {'radius'}, l);
    
    if ~isempty(circleList.item(0).getElementsByTagName('center').item(0))
        x_center = aux_getXEl(circleList, {'center', 'x'}) + x_state;
        y_center = aux_getXEl(circleList, {'center', 'y'}) + y_state;
    else
        x_center = x_state;
        y_center = y_state;
    end 
    
    % convering circle into a polygon
    n = 50; % number of points approximating the circle
    theta = (0:n-1)*(2*pi/n);
    x = x_center + radius*cos(theta);
    y = y_center + radius*sin(theta);
    
    % add all the circles(polygons)
    shape_polygon = union(shape_polygon, polyshape(x,y,'Simplify',false));
end

% polygon
polygonList = aux_getXList(shapeList, {'polygon'});
for l = 0:(polygonList.getLength()-1)
    polygon_pointList = aux_getXList(polygonList, {'point'}, l);
    points = zeros(polygon_pointList.getLength(), 2);
    %y = zeros(1,polygon_pointList.getLength());
    for k = 0:(polygon_pointList.getLength()-1)
        points((k+1), 1) = aux_getXEl(polygon_pointList, {'x'}, k);
        points((k+1), 2) = aux_getXEl(polygon_pointList, {'y'}, k);
    end

    % add all the polyshapes
    shape_polygon = union(shape_polygon, polyshape(points,'Simplify',false));
end

vertices =  shape_polygon.Vertices;
end

function dynamicBuffer = aux_dynamic_fct(obstaclesList, i)
% This function returns a list of dynamic Obstacles extracted from
% obstaclesList.item(i). It handles occupancySets and trajectories.

% initial state
ShapeList_initialstate = aux_getXList(obstaclesList, {'shape'}, i);
stateList_initialstate = aux_getXList(obstaclesList, {'initialState'}, i);
initialstate_vertices = aux_shape_fct(ShapeList_initialstate, stateList_initialstate);

occupancySetXList = aux_getXList(obstaclesList, {'occupancySet'}, i);
trajectoryList = aux_getXList(obstaclesList, {'trajectory'}, i);

% by occupancySet
if ~isempty(occupancySetXList.item(0))
    occupancyList = aux_getXList(occupancySetXList, {'occupancy'});
    
    % creating obstacles
    dynamicBuffer = aux_occupancy_fct(occupancyList);
    
% by trajectory
elseif ~isempty(trajectoryList.item(0))
    trajectorystateList = aux_getXList(trajectoryList, {'state'});
    
    % creating obstacles
    dynamicBuffer = aux_trajectory_fct(trajectorystateList, ShapeList_initialstate);
else
    fprintf(strcat('Found dynamic obstacle without dynamics or with probability',...
        'distribution, which cannot be handled by converter. Resuming...\n'));
end

% inserting the initial state
dynamicBuffer(1).polygon = initialstate_vertices;
dynamicBuffer(1).time = 0;

end

function obsmat = aux_dynamic_mat(obstaclesList, timeStepSize, i)
% This function returns a list of dynamic Obstacles extracted from
% obstaclesList.item(i). It handles occupancySets and trajectories.

% initial state
initialTimeStep = zeros(1,7);
stateList_initialstate = aux_getXList(obstaclesList, {'initialState'}, i);
initialTimeStep(1,1) = i+1;
initialTimeStep(1,2) = aux_getXEl(stateList_initialstate, {'position', 'point', 'x'}); %x_position
initialTimeStep(1,3) = aux_getXEl(stateList_initialstate, {'position', 'point', 'y'}); %y_position
initialTimeStep(1,4) = aux_getXEl(stateList_initialstate, {'orientation', 'exact'}); %orientation
initialTimeStep(1,5) = aux_getXEl(stateList_initialstate, {'velocity', 'exact'}); %velocity

trajectoryList = aux_getXList(obstaclesList, {'trajectory'}, i);
% by trajectory
if ~isempty(trajectoryList.item(0))
    trajectorystateList = aux_getXList(trajectoryList, {'state'});
    trajectorystateListLength = trajectorystateList.getLength;

    obsmat = zeros((trajectorystateListLength + 1),7);
    obsmat(1,:) = initialTimeStep;
    
    if ~isempty(trajectorystateList.item(0))
        for j = 0 : (trajectorystateListLength -1)
            % if mod(j, 6) == 0
                obsmat(j+2,1) = i+1;
    
                midval = aux_getXList(trajectorystateList, {'position', 'point', 'x'}, j); %x_position
                obsmat(j+2,2) = aux_getXEl(midval); %x_position
    
                obsmat(j+2,3) = aux_getXEl(trajectorystateList, {'position', 'point', 'y'}, j); %y_position
                obsmat(j+2,4) = aux_getXEl(trajectorystateList, {'orientation', 'exact'}, j); %orientation
                obsmat(j+2,5) = aux_getXEl(trajectorystateList, {'velocity', 'exact'}, j); %velocity
        
                % Calculate Velocities of all points
                obsmat(j+1,6) = (obsmat(j+2,2) - obsmat(j+1,2)) / timeStepSize; %velocity x
                obsmat(j+1,7) = (obsmat(j+2,3) - obsmat(j+1,3)) / timeStepSize; %velocity y
            % end
        end
    end

    obsmat(:,1) = i+1;
    % Remove last data point as it has velocity = 0
    obsmat(size(obsmat),:) = [];

else
    fprintf(strcat('Found dynamic obstacle without dynamics or with probability',...
        'distribution, which cannot be handled by converter. Resuming...\n'));
end

end

function occupancy = aux_occupancy_fct(occupancyList)
% This function returns a list of occupancy structs for a given xDocList
occupancy(occupancyList.getLength).polygon = [];
for l = 0:(occupancyList.getLength()-1)
    
    % polygon
    shapeList_occupancy = aux_getXList(occupancyList, {'shape'}, l);
    if ~isempty(shapeList_occupancy.item(0))
        occupancy(l+2).polygon = aux_shape_fct(shapeList_occupancy);
    else
        throw(CORAerror('CORA:converterIssue',...
            'No shapes found in occupancy.'));
    end
    
    % time
    timeList_occupancy = aux_getXList(occupancyList, {'time'}, l);
    if ~isempty(timeList_occupancy.item(0).getElementsByTagName('exact').item(0))
        occupancy(l+2).time = aux_getXEl(timeList_occupancy, {'exact'});
        
    elseif ~isempty(timeList_occupancy.item(0).getElementsByTagName('intervalStart').item(0))
        occupancy(l+2).time(1) = aux_getXEl(timeList_occupancy, {'intervalStart'});
        occupancy(l+2).time(2) = aux_getXEl(timeList_occupancy, {'intervalEnd'});
    end
end

end

function trajectory = aux_trajectory_fct(trajectoryStateXList, shapeXList_initialstate)
% This function returns a list of trajectory polygon structs for a given xDocList

trajectory(trajectoryStateXList.getLength).polygon = [];
for l = 0:(trajectoryStateXList.getLength()-1)
    % polygon
    trajectory(l+2).polygon = aux_shape_fct(shapeXList_initialstate, trajectoryStateXList, l);
    
    % time
    timeList = aux_getXList(trajectoryStateXList, {'time', 'exact'}, l);

    if ~isempty(timeList.item(0))
        trajectory(l+2).time = aux_getXEl(timeList);
    end
end
end


function elements = aux_getXEl(xDocList, varargin)
% Returns elements queried in cell array strings. If no vector of indexes
% is provided, the first element (0-element) is returned. Indexes in
% indexes are taken in order, starting with the first element until no
% further elements are provided. There may be more elements in strings than
% in indexes.

[strings,indices] = setDefaultValues({[],0},varargin);
if isempty(indices)
    indices = 0;
end

if ~isempty(strings)
    elements = aux_getXEl(xDocList.item(indices(1)).getElementsByTagName(strings{1}), strings(2:end), indices(2:end));
else
    assert(length(indices) == 1, 'wrong number of indexes provided')
    assert(xDocList.getLength() == 1, 'Several items found, where one was expected.')
    
    elements = str2double(xDocList.item(indices).getTextContent);
end

end

function elements = aux_getXList(xDocList, strings, varargin)
% Returns elements queried in cell array strings. If no vector of indexes
% is provided, the first element (0-element) is returned. Indexes in
% indexes are taken in order, starting with the first element until no
% further elements are provided. There may be more elements in strings than
% in indexes.

indices = setDefaultValues({0},varargin);
if isempty(indices)
    indices = 0;
end

if length(strings) > 1
    elements = aux_getXList(xDocList.item(indices(1)).getElementsByTagName(strings{1}), strings(2:end), indices(2:end));
else
    elements = xDocList.item(indices(1)).getElementsByTagName(strings{1});
end

end

function inter_or_exact = aux_interval_or_exact(XList, varargin)
% This function returns either an interval of two values or a single value
% from item(idx) given in XList, depending on whether the childNodes are
% 'exact' or 'intervalStart'/'intervalEnd'

idx = setDefaultValues({0},varargin);

if ~isempty(XList.item(idx).getElementsByTagName('exact').item(0))
    inter_or_exact = aux_getXEl(XList, {'exact'}, idx);
elseif ~isempty(XList.item(idx).getElementsByTagName('intervalStart').item(0))
    inter_or_exact = interval(aux_getXEl(XList, {'intervalStart'}, idx), aux_getXEl(XList, {'intervalEnd'}, idx));
end

end

% ------------------------------ END OF CODE ------------------------------
