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
    child.Binds(i).values_text = applyRenames(child.Binds(i).values_text,bind.keys,bind.renames);
    % substitute variables for mapped values
    child.Binds(i).values = applySymMapping(child.Binds(i).values,bind.keys,bind.values);
end

end