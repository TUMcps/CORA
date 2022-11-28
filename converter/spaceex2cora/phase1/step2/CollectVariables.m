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

% Author:       ???
% Written:      ???
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% store variables and labels
listOfVars = struct([]);
listOfLabels = struct([]);

% count variables used to index structs above
h_nLab = 0;
h_nVar = 0;

% Iterate over parameters to:
% - split labels and other variables
% - detect local variables

% these attributes are needed to parse a param definition
generalAttr = {'name','type','local'};

% loop over all parameters
for i = 1:length(paramStruct)

    % check existence of parsed fields
    if ~all(isfield(paramStruct{i}.Attributes,generalAttr))
        throw(CORAerror('CORA:converterIssue',...
            ['Parameter ' num2str(i) 'lacking one of the necessary fields: '...
            '"name", "type", "local".']));
    end
    
    param_type = paramStruct{i}.Attributes.type;
    switch param_type

        case 'label'
            % parsing label
            h_nLab = h_nLab +1;
            listOfLabels(h_nLab).name = paramStruct{i}.Attributes.name;

        otherwise
            % 'int','any' currently not recieving special treatment
            if ~strcmp(param_type,'real')
                warning('Parameter %s: type "%s" not supported, treating as "real".',...
                    listOfVars(h_nVar).name,param_type);
            end
            % parsing variable
            h_nVar = h_nVar + 1;
            name_unsafe = paramStruct{i}.Attributes.name;
            
            % Unfortunately, the symbolic Toolbox interprets variables
            % named "i", "j", "I", or "J" as the imaginary number.
            % -> transform all variable names to avoid this.
            listOfVars(h_nVar).name = replaceImagVarnames(name_unsafe);
            
    end
    
    % Other fields are not used currently.
    % Their parsing code was cut, check git history if it's needed again.
end

%------------- END OF CODE --------------
