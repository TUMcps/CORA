function [configParams,configSpecs,spec_mapping] = ...
    parseSpaceExConfig(filename,state_names,input_names,component_names,location_names)
% parseSpaceExConfig - Constructs structs for parameters and specifications
%    for reachability analysis as given in a SpaceEx config-file
%
% Syntax:
%    [configParams,configSpecs,spec_mapping] = ...
%       parseSpaceExConfig(filename,state_names,input_names,component_names,location_names)
%
% Inputs:
%    filename - path to the cfg-file that contains the SpaceEx-configuration 
%    state_names - Struct array containing structs with field "name"
%                  specifying the names of the state variables
%    input_names - Struct array containing structs with field "name"
%                  specifying the names of the input variables
%    component_names - String array containing the names of all automaton
%                      components
%    location_names - Cell array where entry i contains a list of all
%                     location names of automaton component i
%    NOTE: All input parameters containing names must be in order!
%
% Outputs:
%    configParams - Struct object containing all relevant
%                   reachability-analysis parameters given in the 
%                   specified .cfg file
%    configSpecs   - CORA specification object denoting forbidden areas
%                    within the state space, as given in the specified
%                    .cfg file
%    spec_mapping - cell array specifying which specifications are relevant
%                   to which locations of which subcomponent
%                   e.g.: spec_mapping =  {{[1,2]; [1]}, {[2]; [3]}} means
%                     subcomponent 1: specs 1,2 for location 1, spec 1 for
%                     location 2
%                     subcomponent 2: spec 2 for location 1, spec 3 for
%                     location 2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Maximilian Perschl
% Written:       05-September-2021
% Last update:   30-September-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% open file descriptor
configFile = fopen(filename,'r');

% parse arguments

% get first line
currentLine = fgetl(configFile);

while ~feof(configFile)
    % remove spaces
    currentLine = strrep(currentLine,' ','');
    
    % if line is a comment or empty, skip
    if strlength(currentLine) == 0 || currentLine(1) == '#'
        currentLine = fgetl(configFile);
        continue;
    end
    
    % get argument
    splittedLine = strsplit(currentLine,'=');
    
    % continue based on argument: arguments 'initially' may be longer than
    % a single line, therefore its parser function receives the file descriptor
    if strcmp(splittedLine(1),'initially')
        [configParams.R0,configParams.u] = ...
            parseInitial(currentLine,configFile,state_names,input_names);
    elseif strcmp(splittedLine(1),'forbidden')
        [configSpecs, spec_mapping] = ...
            parseForbidden(currentLine,state_names,component_names,location_names);
    elseif strcmp(splittedLine(1),'time-horizon')
        configParams.Tfinal = str2double(splittedLine(2));
    elseif strcmp(splittedLine(1),'sampling-time')
        configParams.Tsample = str2double(splittedLine(2));
    end
    % all other arguments are irrelevant for CORA and can be skipped

    % move on to next line - file pointer is at the beginning of the next relevant line
    currentLine = fgetl(configFile);
end
    

% check if initial set and time horizon given 
hasInitState = isfield(configParams,'R0');
hasTfinal = isfield(configParams,'tFinal');

if ~(hasInitState && hasTfinal)
    disp('The given params-struct is missing crucial arguments!\nPlease add:');
    if ~hasInitState
        disp('  Initial continuous state (R0)');
    end
    if ~hasTfinal
        disp('  Time horizon (tFinal)');
    end
end

fclose(configFile);

% ------------------------------ END OF CODE ------------------------------
