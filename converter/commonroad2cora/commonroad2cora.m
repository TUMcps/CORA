function [statObs,dynObs,x0,goalSet,lanelets,information] = commonroad2cora(filename,varargin)
% commonroad2cora - convert CommonRoad XML-file scenario to CORA model
%
% Syntax:
%    [statObs,dynObs,x0,goalSet,lanelets] = commonroad2cora(filename);
%
% Inputs:
%    filename - name of the CommonRoad XML-file scenario (specified as string)
%    'verbose', true/false - name-value pair specifying whether conversion
%           information should be printed to the terminal. Default is true.
%    'roadBoundary', true/false - name-value pair specifying whether a road
%           boundary should be created as a static obstacle. Default is true.
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
%    [statObs,dynObs,x0,goalSet,lanelets] = ...
%        commonroad2cora('USA_US101-1_1_S-1','verbose',false);
%
%    figure; hold on;
%    for i=1:length(lanelets)
%        plot(lanelets{i}.set,'FaceColor',[.7 .7 .7]);
%    end
%
%    plot(goalSet{1}.set,[1,2],'FaceColor','r','FaceAlpha',0.5);
%    plot(x0.x,x0.y,'.g','MarkerSize',20);
%
%    for i=1:length(dynObs)
%        plot(dynObs{i}.set,[1,2],'FaceColor','b');
%    end
%    for i=1:length(statObs)
%        plot(statObs{i}.set,'FaceColor','y');
%    end
%    % Note: change 'EdgeColor' for static obstacles for better visibility
%    % of road boundary
%
%    xlim([-32,35]); ylim([-14,15]);
%    axis equal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: spaceex2cora

% Authors:       Farah Atour, Niklas Kochdumper, Philipp Gassert
% Written:       28-April-2020
% Last update:   17-September-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

conTimerStart = tic;

% ----------------------------
% ----- Reading options ------
% ----------------------------

defaultVerbose = true;
defaultRoadBoundary = true;

p = inputParser;
addParameter(p, 'verbose', defaultVerbose, @(x) islogical(x) || (x == 1) || (x == 0) );
addParameter(p, 'roadBoundary', defaultRoadBoundary, @(x) islogical(x) || (x == 1) || (x == 0) );

parse(p, varargin{:});

verbose = p.Results.verbose;
roadBoundary = p.Results.roadBoundary;

% ----------------------------
% --- Reading the XML file ---
% ----------------------------

% add '.xml' to the filename if missing
if ~contains(filename, '.xml')
    filename = strcat(filename,'.xml');
end
xDoc = xmlread(filename);

% assert number of scenarios
assert(xDoc.getElementsByTagName('commonRoad').getLength == 1, 'Cannot handle more than one scenario!');

% --- Retreive commonRoad version ---
commonroad_version = xDoc.getElementsByTagName('commonRoad').item(0).getAttribute('commonRoadVersion');

% -----------------------------------------------
% -- Extract all static and dynamic obstacles ---
% -----------------------------------------------

timeStep = str2double(xDoc.item(0).getAttribute('timeStepSize'));
dynamic_counter = 0;        % counts number of dynamic obstacles
static_counter  = 0;        % counts number of static obstacles
total_dynamic_counter = 0;  % counts total number of all (time independet) dynObs
statObs = [];
dynObs = [];
boundaryMessage_flag = 1;

w = warning();
warning('off');

% --- Extracting obstacles for 2020a version ---
if strcmp(commonroad_version, '2020a')
    dynamicObstaclesList = xDoc.getElementsByTagName('dynamicObstacle');
    dynamicObstacles_length = dynamicObstaclesList.getLength;
    staticObstaclesList = xDoc.getElementsByTagName('staticObstacle');
    staticObstacles_length = staticObstaclesList.getLength;
    
    % extracting dynamic obstacles
    for i = 0:(dynamicObstacles_length-1)
        dynamic_counter = dynamic_counter + 1;
        
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
    end
    
    % extracting static obstacles
    for i = 0:(staticObstacles_length-1)
        
        % Inform about existing boundary and possibility to suppress
        if (staticObstaclesList.item(i).getElementsByTagName('type').getLength > 0) &&...
                strcmp(aux_getXEl(staticObstaclesList, {'type'}, i), 'roadBoundary') && roadBoundary && boundaryMessage_flag
            fprintf(strcat('Road Boundary contained in file! Additional Road Boundary is created; ',...
                'to suppress additional creation, call with name-value pair "roadBoundary", false.\n'));
        end
        
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
    
% --- Extracting obstacles from 2018b version ---
elseif strcmp(commonroad_version, '2018b')
    obstaclesList = xDoc.getElementsByTagName('obstacle');
    obstacles_length = obstaclesList.getLength;
    
    if ~isempty(obstaclesList.item(0))
        
        for i = 0:(obstacles_length-1)
                       
            % role (static or dynamic)
            role = char(obstaclesList.item(i).getElementsByTagName('role').item(0).getTextContent);
            
            % dynamic obstacles
            if strcmp(role, 'dynamic')
                dynamic_counter = dynamic_counter + 1;
                
                dynamicBuffer = aux_dynamic_fct(obstaclesList, i);
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
            % static obstacles    
            else
                
                % initial state
                ShapeList_initialstate = aux_getXList(obstaclesList, {'shape'}, i);
                stateList_initialstate = aux_getXList(obstaclesList, {'initialState'}, i);
                initialstate_vertices = aux_shape_fct(ShapeList_initialstate, stateList_initialstate);
                
                % Inform about existing boundary and possibility to suppress
                if (obstaclesList.item(i).getElementsByTagName('type').getLength > 0) &&...
                        strcmp(aux_getXEl(obstaclesList, {'type'}, i), 'roadBoundary') && roadBoundary && boundaryMessage_flag
                    fprintf(strcat('Road Boundary contained in file! Additional Road Boundary is created; ',...
                        'to suppress additional creation, call with name-value pair "roadBoundary", false.\n'));
                end
                
                % split into regions to obtain single-regions obstacles
                initialstate_polyshape = polyshape(initialstate_vertices);               
                reg = regions(initialstate_polyshape);
                
                if length(reg) > 1
                    fprintf('More than one region in static obstacle.\n')
                end
                
                for l = 1:length(reg)
                    static_counter = static_counter + 1;
                    vert = reg(l,1).Vertices;
                    statObs{end+1,1} = polygon(vert(:,1),vert(:,2));
                end
            end
        end
    end
else
    fprintf('Unidentified version (not 2018b or 2020a) when extracting obstacle information. No obstacles extracted.');
end

% -----------------------------------------------
% -- Extract all lanelets of the road network ---
% -----------------------------------------------

laneletList = xDoc.getElementsByTagName('lanelet');

% --- Preallocate for better memory and speed ---
lanelet_length = laneletList.getLength();
lanelet_flag = lanelet_length > 0;
laneletsList = struct('id', cell(lanelet_length,1), ...
    'leftBound', cell(lanelet_length,1), ...
    'rightBound', cell(lanelet_length,1));
lanelets = cell(lanelet_length,logical(lanelet_length));
bound_iterator = 0;

for i = 0:(lanelet_length-1)
    
    % retrieve id and set to -1 for references from goal region specification
    if ~strcmp(laneletList.item(i).getAttribute('id'),'')
        laneletsList(i+1).id = str2double(laneletList.item(i).getAttribute('id'));
    else
        laneletsList(i+1).id = -1;
        continue
    end
    
    % left bound
    leftbound_pointList = aux_getXList(laneletList, {'leftBound', 'point'}, i);
    laneletsList(i+1).leftBound = zeros(2,leftbound_pointList.getLength());
    
    for j = 0:(leftbound_pointList.getLength()-1)
        laneletsList(i+1).leftBound(1,j+1) = aux_getXEl(leftbound_pointList, {'x'}, j);
        laneletsList(i+1).leftBound(2,j+1) = aux_getXEl(leftbound_pointList, {'y'}, j);
    end
    
    % right bound
    rightbound_pointList = aux_getXList(laneletList, {'rightBound', 'point'}, i);
    laneletsList(i+1).rightBound = zeros(2,rightbound_pointList.getLength());
    
    for j = 0:(rightbound_pointList.getLength()-1)
        laneletsList(i+1).rightBound(1,j+1) = aux_getXEl(rightbound_pointList, {'x'}, j);
        laneletsList(i+1).rightBound(2,j+1) = aux_getXEl(rightbound_pointList, {'y'}, j);
    end
    
    % vertices of the left bound lanelet (eases handling below)
    x_left = laneletsList(i+1).leftBound(1,:);
    y_left = laneletsList(i+1).leftBound(2,:);
    
    % vertices of the right bound lanelet (eases handling below)
    x_right = laneletsList(i+1).rightBound(1,:);
    y_right = laneletsList(i+1).rightBound(2,:);
    
    % insertion into output cell array
    lanelets{i+1} = polygon([x_left, flip(x_right)], [y_left, flip(y_right)]);
    if lanelets{i+1}.set.NumRegions > 1
        fprintf(strcat('Warning: Lanelet consists of %d regions!',...
            ' Assuming imprecision contained in pointlist input file.',...
            ' Continuing without alteration.\n'), lanelets{i+1}.set.NumRegions);
    end
    
    % scenarion boundaries for later use in roadBoundary
    if i == 0
        lim_x = [min([x_left, x_right]), max([x_left, x_right])];
        lim_y = [min([y_left, y_right]), max([y_left, y_right])];
    else
        lim_x = [min([lim_x(1), x_left, x_right]), max([lim_x(2), x_left, x_right])];
        lim_y = [min([lim_y(1), y_left, y_right]), max([lim_y(2), y_left, y_right])];
    end
    
    % -----------------------------------------------
    % --- Prepare drivable space for roadBoundary ---
    % -----------------------------------------------
    
    if roadBoundary
        
        % retrieve all left neighbours of lane i
        if laneletList.item(i).getElementsByTagName('adjacentLeft').getLength == 1
            ref_id = str2double(laneletList.item(i).getElementsByTagName('adjacentLeft').item(0).getAttribute('ref'));
            direction = laneletList.item(i).getElementsByTagName('adjacentLeft').item(0).getAttribute('drivingDir');
            
            % find referred lane
            for j = 1:i
                if (laneletsList(j).id == ref_id) && strcmp(direction, 'same')
                    bound = laneletsList(j).rightBound;
                elseif (laneletsList(j).id == ref_id) && strcmp(direction, 'opposite')
                    bound = laneletsList(j).leftBound;
                else
                    continue
                end
                
                bound_iterator = bound_iterator + 1;
                roadBoundaryResiduals(bound_iterator).x = [x_left, bound(1,:)];
                roadBoundaryResiduals(bound_iterator).y = [y_left, bound(2,:)];
                break;
            end
        elseif laneletList.item(i).getElementsByTagName('adjacentLeft').getLength > 1
            throw(CORAerror('CORA:converterIssue',...
                'More than one adjacentLeft lanelet. Violates assumptions; Road Boundary functionality fails.'));
        end
        
        % retrieve all right neighbours of lane i
        if laneletList.item(i).getElementsByTagName('adjacentRight').getLength == 1
            ref_id = str2double(laneletList.item(i).getElementsByTagName('adjacentRight').item(0).getAttribute('ref'));
            direction = laneletList.item(i).getElementsByTagName('adjacentRight').item(0).getAttribute('drivingDir');
            
            % find referred lane
            for j = 1:i 
                if (laneletsList(j).id == ref_id) && strcmp(direction, 'same')
                    bound = laneletsList(j).leftBound;
                elseif  (laneletsList(j).id == ref_id) && strcmp(direction, 'opposite')
                    bound = laneletsList(j).rightBound;
                else
                    continue
                end
                
                bound_iterator = bound_iterator + 1;
                roadBoundaryResiduals(bound_iterator).x = [x_right, bound(1,:)];
                roadBoundaryResiduals(bound_iterator).y = [y_right, bound(2,:)];
                break;
            end
        elseif laneletList.item(i).getElementsByTagName('adjacentRight').getLength > 1
            throw(CORAerror('CORA:converterIssue',...
                'More than one adjacentRight lanelet. Violates assumptions; Road Boundary functionality fails.'));
        end
        
        overshoot_len = 4; % additional free space before and behind road ends (in m)
        
        % retrieve all predecessors of lane i or clear space in front by overshoo_len
        num_pred = laneletList.item(i).getElementsByTagName('predecessor').getLength;
        if ~num_pred
            
            % clear space in front
            bound_iterator = bound_iterator + 1;
            points = [x_left(1), y_left(1); x_right(1), y_right(1)];
            % right angle vector:
            overshoot_vec = ((points(2,:) - points(1,:)) / pdist(points)) * overshoot_len;
            roadBoundaryResiduals(bound_iterator).x = [x_left(1); x_right(1);...
                x_right(1) + overshoot_vec(2); x_left(1) + overshoot_vec(2)];
            roadBoundaryResiduals(bound_iterator).y = [y_left(1); y_right(1);...
                y_right(1) - overshoot_vec(1); y_left(1) - overshoot_vec(1)];
        else
            
            % find referred lane
            for j = 0:(num_pred-1)
                ref_id = str2double(laneletList.item(i).getElementsByTagName('predecessor').item(j).getAttribute('ref'));
                
                for k = 1:i
                    if laneletsList(k).id == ref_id
                        bound_iterator = bound_iterator + 1;
                        roadBoundaryResiduals(bound_iterator).x = [x_left(1); x_right(1);...
                            laneletsList(k).rightBound(1,end); laneletsList(k).leftBound(1,end)];
                        roadBoundaryResiduals(bound_iterator).y = [y_left(1); y_right(1);...
                            laneletsList(k).rightBound(2,end); laneletsList(k).leftBound(2,end)];
                        break;
                    end
                end
            end
        end
        
        % retrieve all successors of lane i or clear space behind by overshoo_len
        num_suc = laneletList.item(i).getElementsByTagName('successor').getLength;
        if ~num_suc
            
            % clear space behind
            bound_iterator = bound_iterator + 1;
            points = [x_left(end), y_left(end);x_right(end), y_right(end)];
            % right angle vector:
            overshoot_vec = ((points(2,:) - points(1,:)) / pdist(points)) * overshoot_len;
            roadBoundaryResiduals(bound_iterator).x = [x_left(end); x_right(end);...
                x_right(end) - overshoot_vec(2); x_left(end) - overshoot_vec(2)];
            roadBoundaryResiduals(bound_iterator).y = [y_left(end); y_right(end);...
                y_right(end) + overshoot_vec(1); y_left(end) + overshoot_vec(1)];
        else
            
            % find referred lanes
            for j = 0:(num_suc-1)
                ref_id = str2double(laneletList.item(i).getElementsByTagName('successor').item(j).getAttribute('ref'));
                
                for k = 1:i
                    if laneletsList(k).id == ref_id
                        bound_iterator = bound_iterator + 1;
                        roadBoundaryResiduals(bound_iterator).x = [x_left(end); x_right(end);...
                            laneletsList(k).rightBound(1,1); laneletsList(k).leftBound(1,1)];
                        roadBoundaryResiduals(bound_iterator).y = [y_left(end); y_right(end);...
                            laneletsList(k).rightBound(2,1); laneletsList(k).leftBound(2,1)];
                        break;
                    end
                end
            end
        end
    end
end

% delete non-lanelet-constituting lanelets (i.e. lanelet references in goal
% regions that were picked up by getElementsByTagName)
for i = lanelet_length:-1:1
    if laneletsList(i).id == -1
        laneletsList(i) = [];
        lanelets(i) = [];
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
            assert(positionList.getLength() < 2, 'Violation of standard: More than one position for goalState.\n')
            if ~isempty(positionList.item(0))
                
                % lanelets
                laneletList = positionList.item(0).getElementsByTagName('lanelet');
                lanelet_length = laneletList.getLength;
                polyshapeBuffer = polyshape();
                for m = 0:(lanelet_length-1)
                    ref = str2double(laneletList.item(m).getAttribute('ref'));
                    idx = [laneletsList.id] == ref;
                    polyshapeBuffer = union(polyshapeBuffer,...
                        polyshape([laneletsList(idx).leftBound, flip(laneletsList(idx).rightBound, 2)]'));
                end
                
                % region by lanelets if found above
                if lanelet_length
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
        fprintf('No planning problem found. Will return empty vectors x0, goalSet');
    end
    x0 = [];
    goalSet = [];
end

% -----------------------------------------------
% -- Extract further information ----------------
% -----------------------------------------------

if strcmp(commonroad_version, '2020a')
    information = struct();
    
    try
        % tags
        tags_all = xDoc.getElementsByTagName('scenarioTags').item(0).getElementsByTagName('*');
        tags = cell(length(tags_all));
        
        for i = 0:length(tags_all)
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

% --- Static Obstacles unit test---
assert(length(statObs) == static_counter, 'Counter mismatch for static obstacles.')

% --- Dynamic Obstacles unit test ---
assert(length(dynObs) == total_dynamic_counter, 'Counter mismatch for dynamic obstacles.')

% --- Road boundary creation ---
if roadBoundary && lanelet_flag
    bounds = [lim_x(1)-1,lim_y(1)-1; lim_x(1)-1,lim_y(2)+1;
              lim_x(2)+1,lim_y(2)+1; lim_x(2)+1,lim_y(1)-1];
    
    roadBoundaryPolygon = polyshape(bounds);
    
    % subtraction of lanelets (twice, as residuals remained for unknown
    % reasons, if only done once)
    for i = [1:length(laneletsList), length(laneletsList):-1:1]
        points = [laneletsList(i).leftBound(1,:), flip(laneletsList(i).rightBound(1,:));...
                laneletsList(i).leftBound(2,:), flip(laneletsList(i).rightBound(2,:))]';
        roadBoundaryPolygon = subtract(roadBoundaryPolygon, polyshape(points));
    end
    
    % subtraction of residuals (twice, as residuals remained for unknown
    % reasons, if only done once)
    for i = [1:bound_iterator, bound_iterator:-1:1]
        roadBoundaryPolygon = subtract(roadBoundaryPolygon,...
            polyshape(roadBoundaryResiduals(i).x, roadBoundaryResiduals(i).y));
    end
    
    buffer = polygon(roadBoundaryPolygon.Vertices(:,1), roadBoundaryPolygon.Vertices(:,2));
    statObs{end+1,1} = buffer;
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
        ' with\n%d lanelets,\n',...
        '%d static obstacles,\n%d dynamic obstacles,\n',...
        '%d planning problem with\n%d goal regions in\n%.2f seconds.\n--------------\n'),...
        length(lanelets), static_counter, dynamic_counter, egoVehiclesList.getLength, num_goalRegions, converterTime);
end

end


% Auxiliary functions -----------------------------------------------------

function vertices = aux_shape_fct(shapeList, stateXList, idx)
% This function returns vertices for the combined polyshape which is the
% union of all shapes contained in shapeList, possibly combined with
% position and orientation from the optional parameter stateXList. You can
% choose which states from stateXList with index idx are relevant by
% passing the optional parameter idx.

if ~exist('idx', 'var')
    idx = 0;
end

if ~exist('stateXList', 'var')
    orientation_state = 0;
    x_state = 0;
    y_state = 0;
else
    assert(stateXList.getLength() >= idx + 1,...
        sprintf('Expected state list of min length %i; got length %i.\n', idx, stateXList.getLength()));
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
% This function returns a list of dynamic Obstalces extracted from
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

function occupancy = aux_occupancy_fct(occupancyList)
% This function returns a list of occupancy structs for a given xDocList
occupancy(occupancyList.getLength).polygon = [];
for l = 0:(occupancyList.getLength()-1)
    
    % polygon
    shapeList_occupancy = aux_getXList(occupancyList, {'shape'}, l);
    if ~isempty(shapeList_occupancy.item(0))
        occupancy(l+2).polygon = aux_shape_fct(shapeList_occupancy);
    else
        assert(0, 'This is a temporary developer warning indicating unexpected behaviour: No shapes found in occupancy.')
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

function elements = aux_getXEl(xDocList, strings, indexes)
% Returns elements queried in cell array strings. If no vector of indexes
% is provided, the first element (0-element) is returned. Indexes in
% indexes are taken in order, starting with the first element until no
% further elements are provided. There may be more elements in strings than
% in indexes.

if ~exist('strings', 'var')
    strings = [];
end

if ~exist('indexes', 'var') || isempty(indexes)
    indexes = 0;
end

if ~isempty(strings)
    elements = aux_getXEl(xDocList.item(indexes(1)).getElementsByTagName(strings{1}), strings(2:end), indexes(2:end));
else
    assert(length(indexes) == 1, 'wrong number of indexes provided')
    assert(xDocList.getLength() == 1, 'Several items found, where one was expected.')
    
    elements = str2double(xDocList.item(indexes).getTextContent);
end

end

function elements = aux_getXList(xDocList, strings, indexes)
% Returns elements queried in cell array strings. If no vector of indexes
% is provided, the first element (0-element) is returned. Indexes in
% indexes are taken in order, starting with the first element until no
% further elements are provided. There may be more elements in strings than
% in indexes.

if ~exist('indexes', 'var') || isempty(indexes)
    indexes = 0;
end

if length(strings) > 1
    elements = aux_getXList(xDocList.item(indexes(1)).getElementsByTagName(strings{1}), strings(2:end), indexes(2:end));
else
    elements = xDocList.item(indexes(1)).getElementsByTagName(strings{1});
end

end

function inter_or_exact = aux_interval_or_exact(XList, idx)
% This function returns either an interval of two values or a single value
% from item(idx) given in XList, depending on whether the childNodes are
% 'exact' or 'intervalStart'/'intervalEnd'

if ~exist('idx', 'var')
    idx = 0;
end

if ~isempty(XList.item(idx).getElementsByTagName('exact').item(0))
    inter_or_exact = aux_getXEl(XList, {'exact'}, idx);
elseif ~isempty(XList.item(idx).getElementsByTagName('intervalStart').item(0))
    inter_or_exact = interval(aux_getXEl(XList, {'intervalStart'}, idx), aux_getXEl(XList, {'intervalEnd'}, idx));
end

end

% ------------------------------ END OF CODE ------------------------------
