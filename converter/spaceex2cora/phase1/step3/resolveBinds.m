function instances = resolveBinds(instances)
% resolveBinds - removes constants from listOfVar and writes them into a
%    new field called listOfConstants; this procedure is important to
%    avoid identifying constants as states later on
%
% Syntax:
%    instances = resolveBinds(instances)
%
% Inputs:
%    instances - instantiated network and base components
%
% Outputs:
%    instances - instantiated network and base components, where the latter
%                correctly separate variables and constants
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       13-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of instances
nrInstances = length(instances);

% loop over all instances
for i=1:nrInstances

    % skip all network components
    if ~instances{i}.isNetwork

        % initialize list of constants with empty list
        instances{i}.listOfConstants = struct('name',cell(1,0),'value',cell(1,0));

        % get path of parents of base component
        parentPath = aux_parentPath(instances,i);

        % follow path of parents to check if any of the entries in
        % listOfVar end up being constants
        for j=1:length(parentPath)
            % current parent
            curr_parent = instances{parentPath(j)};
            % current child
            if j == 1
                curr_child = i;
            else
                curr_child = parentPath(j-1);
            end

            % index of current child in children of current parent
            % note: take first child (multiple children aka base components
            % in one network component have same binds)
            idx_child = find(curr_parent.children == curr_child,1,'first');

            % read out bind from parent to current child
            bind = curr_parent.Binds(idx_child);

            % list of variables in base component
            idx_keep = true(numel(instances{i}.listOfVar),1);

            % loop over list of variables, resolve all binds
            for k=1:numel(instances{i}.listOfVar)
                renameIdx = bind.keys == instances{i}.listOfVar(k).name;
            
                if any(renameIdx)
                    % rename state/constant
                    instances{i}.listOfVar(k).name = bind.renames(renameIdx);
            
                    % check if variable is mapped to a constant (converting
                    % a string that does not represent a double value to a
                    % double results in NaN)
                    temp = str2double(bind.renames(renameIdx));
                    if ~isnan(temp)
                        % add to list of constants
                        instances{i}.listOfConstants = [instances{i}.listOfConstants, ...
                            struct('name',bind.keys{renameIdx},'value',temp)];
                        % remove from list of variables
                        idx_keep(k) = false;
                    end
                end

            end

            % re-assign list of variables
            instances{i}.listOfVar = instances{i}.listOfVar(:,idx_keep);

        end

    end

end

end


% Auxiliary functions -----------------------------------------------------

function parentPath = aux_parentPath(instances,idx)
% find path (ordered list of indices) to higher-level network components
% from base component with index idx in cell-array of instances;
% caution: done recursively!

% init parent path
parentPath = [];

for i=1:length(instances)
    if isfield(instances{i},'children') ...
            && any(instances{i}.children == idx)
        parentPath = [parentPath; i; aux_parentPath(instances,i)];
        break
    end
end

end

% ------------------------------ END OF CODE ------------------------------
