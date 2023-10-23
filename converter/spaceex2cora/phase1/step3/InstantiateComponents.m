function instances = InstantiateComponents(templates,rootIdx)
% InstantiateComponents - builds a tree of component instances, beginning
%    with a root component. Instances are created, by applying variable
%    mappings on their templates. Instantiate through deeper trees with
%    breadth-first search, by keeping track of all still open branches.
%
% Syntax:
%    instances = InstantiateComponents(templates,rootIdx)
%
% Inputs:
%    templates - struct containing all base/network components and their
%                properties (variable names, binds, flow equations,
%                invariants, guard sets, reset function, synchronization
%                labels, etc.), as processed from SpaceEx xml-file in the
%                function parseTemplates
%    rootIdx - ID of the root component (default: top-most network component)
%
% Outputs:
%    instances - instantiated components
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

% keep track of instance count and tree depth
num_instances = 1;
depth = 0;

% 1st BFS iteration: create root node
% (the root instance does not need to be altered)
instances = templates(rootIdx);
% init naming system by setting root name
instances{1}.name = instances{1}.id;

% Keep track of all non-leaves/network-components, that were instantiated
% in the last BFS-iteration.
if instances{1}.isNetwork
    branches = 1;
else
    branches = [];
end
% Their children will be instantiated in the next iteration. 

% loop over elements until the entire tree has been explored and all
% components have been instantiated
while ~isempty(branches)
    % increment tree depth
    depth = depth + 1;
    % keep track of open branches created in this iteration
    newBranches = [];
    
    % expand open branches
    for br = 1:length(branches)
        
        % store the parent instance for quick access
        branch = instances{branches(br)};
        % keep track of created children
        children = [];
        
        % loop over all instantiated components (aka binds)
        for bi = 1:length(branch.Binds)
            % create a child instance for every instantiated component
            
            % update newest-component index
            num_instances = num_instances + 1;
            
            % store bind and template structs for quick access
            bind = branch.Binds(bi);
            childTemplate = templates{bind.idx};
            
            if childTemplate.isNetwork
                % instantiate a network component
                instances{num_instances} = InstantiateNC(branch,childTemplate,bind);
                % remember currently open branch for next iteration
                newBranches = [newBranches,num_instances];
            else
                % instantiate a base component
                instances{num_instances} = InstantiateBC(branch,childTemplate,bind);
            end
            
            % store child references for parent
            children = [children,num_instances];
        end
    
        % add child references to parent
        instances{branches(br)}.children = children;
    end
    
    % prepare the next layer of branches for instantiation
    branches = newBranches;

end
% if "branches" is empty, no more nodes need to be expanded
% -> BFS traversal complete

% resolve binds in listOfVar: since constants can be inherited over
% multiple levels, we need to track every variable's bind from the base
% level to the highest network component to find out if any bind along the
% way leads to a constant; the heritage can be tracked using the field
% instances.children
instances = resolveBinds(instances);


disp("traversal complete!");
fprintf("  tree depth: %i\n",depth);
fprintf("  total instances: %i\n",num_instances);

% ------------------------------ END OF CODE ------------------------------
