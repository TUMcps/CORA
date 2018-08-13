function instances = InstantiateComponents(templates,rootIdx)
% Build a tree of component instances, beginning with a root component.
% Instances are created, by applying variable mappings on their templates.
% Instanciate through deeper trees with breadth-first search (BFS),
% by keeping track of all still open branches.

% keep track of instance count and tree depth
num_instances = 1;
depth = 0;

% 1st BFS iteration: create root node
% (the root instance does not need to be altered)
instances = {templates{rootIdx}};
% init naming system by setting root name
instances{1}.name = instances{1}.id;

% Keep track of all non-leaves/network-components, that were instanciated
% in the last BFS-iteration.
if instances{1}.isNetwork
    branches = [1];
else
    branches = [];
end
% Their children will be instanciated in the next iteration. 

while(~isempty(branches))
    % increment tree depth
    depth = depth + 1;
    % keep track of open branches created in this iteration
    newBranches = [];
    
    % expand open branches
    for br = 1:length(branches)
        
        % store the parent instance for quick access
        branch = instances{br};
        % keep track of created children
        children = [];
        
        for bi = 1:length(branch.Binds)
            % create a child instance for every Bind
            
            % update newest-component index
            num_instances = num_instances + 1;
            
            % store bind and template structs for quick access
            bind = branch.Binds(bi);
            childTemplate = templates{bind.idx};
            
            if childTemplate.isNetwork
                % instantiate a network component
                instances{num_instances} = InstantiateNC(branch,childTemplate,bind);
                % remember open branch for next iteration
                newBranches = [newBranches,num_instances];
            else
                % instantiate a base component
                instances{num_instances} = InstantiateBC(branch,childTemplate,bind);
            end
            
            % store child references for parent
            children = [children,num_instances];
        end
    
        % add child references to parent
        instances{br}.children = children;
    end
    
    % prepare the next layer of branches for instanciation
    branches = newBranches;
end
% if "branches" is empty, no more nodes need to be expanded
% BFS traversal complete

disp("traversal complete")
fprintf("tree depth: %i\n",depth);
fprintf("total instances: %i\n",num_instances);

end

