function [listOfVar, listOfLab ] = CollectVariables(paramStruct)

%INPUT:
%   paramStruct: parameters in SX formate given in a matlab struct
%OUTPUT:
%   listOfVar: list of variables
%   listOfLab: list of labels

% store variables and labels
listOfVar = struct([]);
listOfLab = struct([]);

% count variables used to index structs above
h_nLab = 0;
h_nVar = 0;


% Iterate over parameters to:
%   -split labels and other variables
%   -

% these attributes are needed to parse a param definition
generalAttr = {'name','type'};

for i = 1:length(paramStruct)
    %checking existance of parsed fields
    attrIdx = isfield(paramStruct{i}.Attributes,generalAttr);
    if ~all(attrIdx)
        error('parameter %d lacking one of necessary fields "name","type"',i);
    end
    
    param_type = paramStruct{i}.Attributes.type;
    switch param_type
        case 'label'
            % parsing label
            h_nLab = h_nLab +1;
            listOfLab(h_nLab).name = paramStruct{i}.Attributes.name;
        otherwise
            % 'int','any' currently not recieving special treatment
            if ~strcmp(param_type,'real')
                warning('Parameter %s: type "%s" not supported, treating as "real".'...
                    ,listOfVar(h_nVar).name,param_type);
            end
            % parsing variable
            h_nVar = h_nVar +1;
            name_unsafe = paramStruct{i}.Attributes.name;
            
            % Unfortunately, the symbolic Toolbox interprets variables named "i","j",
            % "I" or "J" as the imaginary number.
            % We perform a transformation on all variable names to avoid this.
            listOfVar(h_nVar).name = replaceImagVarnames(name_unsafe);
    end
    
    % Other fields are not used currently.
    % Their parsing code was cut, check git history if it's needed again.
end

end
