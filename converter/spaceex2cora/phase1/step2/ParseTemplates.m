function [componentTemplates,templateIDs] = ParseTemplates(sxStruct)
% ParseTemplates - parses the struct generated from reading the SpaceEx 
%    xml-file and returns a list of component templates in StructHA format.
%    Additionally, a list of template names is returned (for quick lookup).
%    The order of templates reflects the order in the xml-file. References
%    to other components are stored as indices in the resulting array.
%
% Syntax:
%    [componentTemplates,templateIDs] = ParseTemplates(sxStruct)
%
% Inputs:
%    sxStruct - struct containing all information contained in the SpaceEx
%               xml-file
%
% Outputs:
%    componentTemplates (struct array) - defining the components of the
%       model. Equations/Assignments etc are specified by symbolic objects,
%       variables which are named in the component but defined somewhere
%       else (inputs/constants etc.) are not yet resolved; with fields
%           .id = name of base/network component
%           .isNetwork = true if NC, false if BC
%           .listOfVar.name = name of all variables in the component
%           .States.id/.name (only BC) = meta-data
%           .States.Trans (only BC) = outgoing transitions
%           .States.Flow (only BC) = flow equation
%           .States.Invariant (only BC) = invariant
%           .Binds (only NC) = binds of variables to base components
%    templateIDs (string array) - names of the components in SpaceEx model
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% remove redundant depths of input argument
sxStruct = sxStruct.sspaceex{1};

% number of base + network components in SpaceEx model
num_templates = length(sxStruct.component);
% pre-allocate cell array for templates
componentTemplates = cell(1,num_templates);

% store component names for fast resolution of references
templateIDs = strings(1,num_templates);

% loop over all component definitions (-> templates) in the xml-file,
% process the contained data to convenient formats and store the resulting
% structs in "componentTemplates"
for i = 1:num_templates

    % read i-th component
    currentComp = sxStruct.component{i};

    % components can be "base" components, containing dynamics
    % they are identified by uniquely containing the field '.location'
    if isfield(currentComp,'location')
        fprintf("parsing BaseComponent '%s'...\n",currentComp.Attributes.id);
        parsedComp = ComputeBaseComponent(currentComp);
        parsedComp.isNetwork = false;

    % otherwise, components are "network" components, containing binds,
    % i.e., instantiation of base component templates 
    % they are identified by uniquely containing the field '.bind'
    elseif isfield(currentComp,'bind')
        fprintf("parsing NetworkComponent '%s'...\n",currentComp.Attributes.id);
        % represent references to other templates by their index in this array
        parsedComp = ComputeNetworkComponent(currentComp,templateIDs(1:i-1));
        parsedComp.isNetwork = true;
        
    else
        throw(CORAerror('CORA:converterIssue',...
            ['Component ' num2str(i) ' seems to be neither base nor network.']));
    end
    
    % store in output arguments
    componentTemplates{i} = parsedComp;
    templateIDs(i) = parsedComp.id;
end

fprintf("done parsing: %i components parsed\n",length(templateIDs));

% ------------------------------ END OF CODE ------------------------------
