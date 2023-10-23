function [listOfVars,listOfLabels] = CollectVariables(paramStruct)
% CollectVariables - Collect list of variable names and synchronization
%    labels occurring in a SpaceEx component (both base and network)
%
% Syntax:
%    [listOfVars,listOfLabels] = CollectVariables(paramStruct)
%
% Inputs:
%    paramStruct (struct) - parameters in SpaceEx-format
%
% Outputs:
%    listOfVars (struct) - list of variables
%    listOfLabels (struct) - list of labels
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

% init struct for variables and labels
listOfVars = struct('name',cell(1,0));
listOfLabels = struct('name',cell(1,0));

% count variables used to index structs above
nrLabels = 0;
nrVars = 0;

% Iterate over parameters to:
% - split labels and other variables
% - detect local variables

% loop over all parameters
for i = 1:length(paramStruct)

    % check existence of parsed fields, following attributes are needed:
    % 'name', 'type', 'local'
    if ~all(isfield(paramStruct{i}.Attributes,{'name','type','local'}))
        throw(CORAerror('CORA:converterIssue',...
            ['Parameter ' num2str(i) 'lacking one of the necessary fields: '...
            '"name", "type", "local".']));
    end
    
    % read out type of parameter
    param_type = paramStruct{i}.Attributes.type;

    % different handling for labels and variables/constants
    switch param_type

        case 'label'
            % parsing synchronization label
            nrLabels = nrLabels + 1;
            listOfLabels(nrLabels).name = paramStruct{i}.Attributes.name;

        otherwise
            % 'int','any' currently not recieving special treatment
            if ~strcmp(param_type,'real')
                warning('Parameter %s: type "%s" not supported, treating as "real".',...
                    listOfVars(nrVars).name,param_type);
            end
            % parsing variable
            nrVars = nrVars + 1;
            name_unsafe = paramStruct{i}.Attributes.name;
            
            % Unfortunately, the symbolic Toolbox interprets variables
            % named "i", "j", "I", or "J" as the imaginary number.
            % -> transform all variable names to avoid this.
            listOfVars(nrVars).name = replaceImagVarnames(name_unsafe);
            
    end
end

% ------------------------------ END OF CODE ------------------------------
