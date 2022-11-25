function comp = AssignTemplate( template, bind, listOfVar, listOfLab)

% INPUT:
%     template: the component template being bound
%     map: struct specifying variable mappings
%     listOfVar: list of Variables of the binding component
%     listOfLab: list of Labels for the template (unused)
% OUTPUT:
%     comp: the component in structHA format with assigned Variables and
%           Constants

comp = template;

if isfield(bind,'map')
    map = bind.map;
else
    map = [];
end
% iterate over variable mappings
for i = 1:length(map)
    
    % 'var' is the variable we want to bind
    var = map{i}.Attributes.key;
    % 'value' is the value we want to assign to 'var'
    value = map{i}.Text;
    
    isLabel = false;
    %Search in the list of Labels
    for j = 1:length(comp.h_listOfLab)
        if strcmp(comp.h_listOfLab(j).name, var)
            % store value in the label
            comp.h_listOfLab(j).value = value;
            isLabel = true;
        end
    end
    if isLabel
        % mapping for labels is finished here
        continue;
    end
    
    % Search for the template variable to be assigned
    varIdx = 0;
    for j = 1:length(comp.h_listOfVar)
        if strcmp(comp.h_listOfVar(j).name, var)
            % remember index
            varIdx = j;
        end
    end
    if varIdx == 0
        error('Error while mapping variables: key "%s" not found.',var);
    end
    
    % Test wether a numeric value is mapped
    if ~isempty(str2num(value))
        comp.h_listOfVar(varIdx).value = value;
        % bound constant mapping is finished here
        continue;
    end
    
    % Search variable being assigned to var
    valueIdx = 0;
    for j = 1:length(listOfVar)
        if strcmp(listOfVar(j).name, value)
            % remember index
            valueIdx = j;
        end
    end
    if valueIdx == 0
        error('Error while mapping variables: assigned variable "%s" not found.',value);
    end
    
    % store name of mapped variable
    listOfVar(valueIdx).mappedTo = value;
    
    % mapping a constant to something makes it a constant
    if listOfVar(valueIdx).constant
        comp.h_listOfVar(varIdx).constant = true;
    end
    
end

end

