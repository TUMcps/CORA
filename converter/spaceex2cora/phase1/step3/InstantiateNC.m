function [child] = InstantiateNC(parent,template,bind)
% Apply a bind to a Network-Component template
% INPUTS (all structs):
%   parent:     instantiated parent component
%   template:   template of BC to instantiate
%   bind:       bind to be used for instantiation

child = template;

% save global name for new instance
child.name = parent.name + "." + bind.localName;

fprintf("instantiating Network Component: %s\n",child.name);

numBinds = length(child.Binds);
for i = 1:numBinds
    % if values_text is a variable name, apply its rename
    %child.Binds(i).values_text = applyRenames(child.Binds(i).values_text,bind.keys,bind.renames);
    % substitute variables for mapped values
    child.Binds(i).values = applySymMapping(child.Binds(i).values,bind.keys,bind.values);
end

for v=1:numel(child.listOfVar)
    var = child.listOfVar(v).name;
    if isempty(bind.renames(bind.keys==var))
        error("Rename error in spaceex2cora-InstantiateNC.m!");
    end
    child.listOfVar(v).name = bind.renames(bind.keys == var);
    for b=1:numel(child.Binds)
        for r=1:numel(child.Binds(b).renames)
            if child.Binds(b).renames(r) == var
               child.Binds(b).renames(r) = bind.renames(bind.keys == var);
            end
        end
    end
end


end