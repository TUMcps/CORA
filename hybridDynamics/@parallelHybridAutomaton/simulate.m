function [obj] = simulate(obj,options)
% simulate - simulates a hybrid automaton
%
% Syntax:  
%    [obj] = simulate(obj, options)
%
% Inputs:
%    obj - hybrid automaton object
%    options - simulation options
%
% Outputs:
%    obj - hybrid automaton object with stored simulation results
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Victor Charlent, Johann Schöpfer, Niklas Kochdumper
% Written: 24-May-2016  
% Last update: 04-July-2018
% Last revision: ---

%------------- BEGIN CODE --------------


% load data from options
tFinal = options.tFinal;                    % final time
finalLoc = cell2mat(options.finalLoc);      % final location

% initialize variables
tInter = options.tStart;    % intermediate time at transitions
loc = options.startLoc;     % actual location
locMat = cell2mat(loc);     % actual location as vector
xInter = options.x0;        % intermediate state at transitions
t = [];         % time vector
x = [];         % state vector
count = 1;      % transition counter


% determine which inputs are specified locally and which globally
inputMap = parseInputMap(obj,options);


% loop over the different locations 
while (tInter<tFinal) && (~isempty(loc)) && ~all(finalLoc == locMat) ... 
        || (count==1)

    % construct new location with local Automaton Product
    currentLocation = locationProduct(obj,loc);

    % construct system input for this location
    options.u = mergeInputVector(obj,inputMap,loc,options.uGlob,options.uCompLoc);
    
    % simulate within the actual location
    [tNew,xNew,locNew,xInter] = simulate(currentLocation,options,tInter,tFinal,xInter);

    % new intermediate time is last simulation time
    tInter = tNew(end);
    
    % store results
    t{count} = tNew;
    x{count} = xNew;
    locList{count} = loc;
    
    % increase counter for transitions
    count = count + 1;
    loc = locNew;
    locMat = cell2mat(loc);
end


% save results to object structure
if isempty(obj.result) || ~isfield(obj.result,'simulation')
    obj.result.simulation{1}.t = t;
    obj.result.simulation{1}.x = x;
    obj.result.simulation{1}.location = locList;
else
    obj.result.simulation{end+1}.t = t;
    obj.result.simulation{end}.x = x;
    obj.result.simulation{end}.location = locList;
end


%------------- END OF CODE --------------

