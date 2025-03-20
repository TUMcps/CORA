function cora2spaceex(Obj,filename,filepath)
% cora2spaceex - convert CORA models to SpaceEx models
%
% Syntax:
%    cora2spaceex(Obj,filename);
%
% Inputs:
%    Obj       - Generated Object
%    filename  - name of the generated SpaceEx model (specified as string)
%    filepath  - desired directory to save the xml-file (optional)
%
% Outputs:
%    xmlFile   - xml-file that contains the SpaceEx-model
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   08-June-2020 (MW, add filepath)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 3
    % if no resultpath is given, use default: "cora/models/CoraConverted"
    filepath = [CORAROOT filesep 'models' filesep 'CoraConverted' filesep];
end

% Headline of the xml-file
% Create the document node(docNode) and the first child element(component)
[docNode, component] = priv_headline_cora2spaceex;

if isa(Obj,'hybridAutomaton')
    
    docNode = priv_hybrid_cora2spaceex(Obj, component, docNode);
    
elseif isa(Obj,'linearSys') || isa(Obj,'nonlinearSys')
    
    % First child element (parameter) of the element node (component)
    priv_parameter_cora2spaceex(Obj, component, docNode);
    % Second child element (location) of the element node (component)
    location = priv_location_cora2spaceex(Obj,component, docNode);
    % First child element (invariant) of the element node (location)
    priv_invariant_cora2spaceex(Obj,location, docNode);
    % Second child element (flow) of the element node (location)
    priv_flow_cora2spaceex(Obj, location, docNode);
    
else
    throw(CORAerror('CORA:notSupported','Given contDynamics class not supported.'));
end

% make directory if it does not exist
if ~exist(filepath,'dir')
    mkdir(filepath);
end
% add to MATLAB path
addpath(filepath);

% concatenate full filename
fname = strcat(filepath,'/',filename,'.xml');
% Writing XML File
xmlwrite(fname,docNode);

end

% ------------------------------ END OF CODE ------------------------------
