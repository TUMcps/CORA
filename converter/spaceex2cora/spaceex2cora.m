function spaceex2cora(xmlFile,varargin)
% spaceex2cora - converts SpaceEx models to CORA models 
%
% Syntax:  
%    spaceex2cora(xmlFile)
%    spaceex2cora(xmlFile,save2pha,rootID,outputName,outputDir)
%    spaceex2cora(xmlFile,save2pha,rootID,outputName,outputDir,cfgFile)
%
% Inputs:
%    xmlFile - path to the xml-file that contains the SpaceEx-model
%    save2pha - conversion to parallel hybrid automaton (default: true)
%               conversion to flat hybrid automaton (false)
%    rootID -  ID of SpaceEx component to be used as root component
%              (specified as a string)
%    outputName - name of the generated CORA model (specified as string)
%    outputDir - path to the desired output directory where all generated
%                files are stored
%    cfgFile - path to the cfg-file that contains the SpaceEx-configuration
%
% Outputs:
%   ---
%
% Example: 
%    spaceex2cora('build_48.xml','sys');
%    cd([coraroot filesep 'models' filesep 'SpaceExConverted' filesep 'build_48']);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      03-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% decision whether parallel HA or flat HA
save2pha = true; % default: parallel HA
if nargin >= 2
    save2pha = varargin{1};
end

% parse the SpaceEx-model
if nargin >= 4
    data = SX2structHA(xmlFile,save2pha,varargin{2},varargin{3});
elseif nargin >= 3
    data = SX2structHA(xmlFile,save2pha,varargin{2});
else
    data = SX2structHA(xmlFile,save2pha);
end

% generate the CORA model files
if nargin >= 4 && ~isempty(varargin{3})
   if nargin >= 5 && ~isempty(varargin{4})
      StructHA2file(data,varargin{3},varargin{4});
   else
      StructHA2file(data,varargin{3});
   end
else
   % default name == xml-file name
   [~,name,~] = fileparts(xmlFile); 
   
   % generate files
   StructHA2file(data,name);
end

% Check if a configuration file was given
if nargin <= 5
    return
end

% read configuration file
cfgFile = varargin{5};

% query model data for information relevant to configuration parsing
% required data is names of components, locations, inputs and states
state_names = [];
input_names = [];
component_names = [];
location_names = {};
for i = 1:length(data.Components)
    % Get names of states, names and components
    state_names = [state_names data.Components{i}.states];
    % exclude dummy inputs
    if ~contains(data.Components{i}.inputs.name,"Dummy")
        input_names = [input_names data.Components{i}.inputs];
    end
    component_names = [component_names data.Components{i}.name];
    % Create empty list of locations for every component
    location_names{i} = [];
    % Fill list of locations
    for j = 1:length(data.Components{i}.States)
        location_names{i} = [location_names{i} data.Components{i}.States(j).name];
    end
end

% parse spaceEx configuration file
[configParams,configSpecs,spec_mapping] = ...
    parseSpaceExConfig(cfgFile,state_names,input_names,component_names,location_names);

% if the converted automaton is flat, we change the format of the
% specification mapping to fit the reach function of flat hybrid
% automata
if ~save2pha
    spec_mapping = spec_mapping{1};
end

% save parameters and specifications
% path is "cora/models/SpaceExConverted/configurationResults"
if ~isempty(varargin{4}) && strlen(varargin{4}) ~= 0
    configurationPath = varargin{4};
    if strlen(varargin{3}) ~= 0
        configurationPath = strcat(configurationPath,varargin{3},"_config.mat");
    else
        configurationPath = strcat(configurationPath,name,"_config.mat");
    end
else
    configurationPath = [CORAROOT filesep 'models' filesep 'SpaceExConverted' filesep name '_config.mat'];
end
save(configurationPath,'configParams','configSpecs','spec_mapping');
    
%------------- END OF CODE --------------