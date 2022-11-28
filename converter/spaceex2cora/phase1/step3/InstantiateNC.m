function child = InstantiateNC(parent,template,bind)
% InstantiateNC - Apply a bind to a Network-Component template
%
% Syntax:
%    child = InstantiateNC(parent,template,bind)
%
% Inputs:
%    parent (struct) - instantiated parent component
%    template (struct) - template of NC to instantiate
%    bind (struct) - bind to be used for instantiation
%
% Outputs:
%    child (struct) - instantiated network component (child of "parent"
%                     network component)
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

% instantiate child from template
child = template;

% save global name for new instance
child.name = parent.name + "." + bind.localName;

% output on command window
fprintf("instantiating Network Component: '%s'\n",child.name);

numBinds = length(child.Binds);
for i=1:numBinds
    % if values_text is a variable name, apply its rename
%     child.Binds(i).values_text = applyRenames(child.Binds(i).values_text,...
%         bind.keys,bind.renames);
    % substitute variables for mapped values
    child.Binds(i).values = applySymMapping(child.Binds(i).values,...
        bind.keys,bind.values);
end

% loop over all variables in network component
for v=1:numel(child.listOfVar)

    % read out variable name
    var = child.listOfVar(v).name;

    % assert that there is a bind for this variable name
    if isempty(bind.renames(bind.keys==var))
        throw(CORAerror('CORA:converterIssue','Renaming error of variables'));
    end

    % rename variable name according to bind
    var_rename = bind.renames(bind.keys == var);
    child.listOfVar(v).name = var_rename;

    % loop over number of instantiated components in given network
    % component; overwrite renamed variable in bind
    for b=1:numBinds
        % index of variable to be renamed
        renameIdx = child.Binds(b).renames == var;
        if ~isempty(renameIdx)
            child.Binds(b).renames(renameIdx) = var_rename;
        end
    end

end

%------------- END OF CODE --------------
