function obj = readONNXNetwork(file_path, varargin)
% readONNXNetwork - reads and converts a network saved in onnx format
%    Note: If the onnx network contains a custom layer, this function will
%    create a +CustomLayer package folder containing all custom layers in
%    your current MATLAB directory.
%
% Syntax:
%    res = neuralNetwork.readONNXNetwork(file_path)
%
% Inputs:
%    file_path - path to file
%    verbose - bool if information should be displayed
%    inputDataFormats - dimensions of input e.g. 'BC' or 'BSSC'
%    outputDataFormats - see inputDataFormats
%    targetNetwork - ...
%    containsCompositeLayers - there are residual connections in the
%    network
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner, Lukas Koller
% Written:       30-March-2022
% Last update:   07-June-2022 (specify in- & outputDataFormats)
%                30-November-2022 (removed neuralNetworkOld)
%                13-February-2023 (simplified function)
%                21-October-2024 (clean up DLT function call)
%                25-June-2025 (LK, use digraph for parsing composite layers)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate parameters
narginchk(1,6)
[verbose, inputDataFormats, outputDataFormats, targetNetwork, ...
    containsCompositeLayers] = setDefaultValues({false, 'BC', 'BC', ...
        'dagnetwork', false}, varargin);

% valid in-/outputDataFormats for importONNXNetwork
validDataFormats = {'','BC','BCSS','BSSC','CSS','SSC','BCSSS','BSSSC', ...
    'CSSS','SSSC','TBC','BCT','BTC','1BC','T1BC','TBCSS','TBCSSS'};
inputArgsCheck({ ...
    {verbose, 'att', 'logical'}; ...
    {inputDataFormats, 'str', validDataFormats}; ...
    {outputDataFormats, 'str', validDataFormats}; ...
    {targetNetwork, 'str', {'dagnetwork', 'dlnetwork'}}; ...
})


if verbose
    disp("Reading network...")
end

% try to read ONNX network using dltoolbox
try
    dlt_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork);

catch ME
    if strcmp(ME.identifier, 'MATLAB:javachk:thisFeatureNotAvailable') && ...
            contains(ME.message,'Swing is not currently available.')
        % matlab tries to indent the code of the generated files for 
        % custom layers, for which (somehow?) a gui is required.
        % As e.g. docker runs don't have a gui, we just try to remove the
        % 'indentcode' function call here ...
        aux_removeIndentCodeLines(ME);

        % re-read network
        dlt_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork);

    else
        rethrow(ME)
    end
end

if containsCompositeLayers
    % Combine multiple layers into blocks to realize residual connections and
    % parallel computing paths.
    layers = aux_groupCompositeLayers(dlt_net.Layers,dlt_net.Connections);
else
    layers = num2cell(dlt_net.Layers);
end

% convert DLT network to CORA network
% obj = neuralNetwork.convertDLToolboxNetwork(dlt_net.Layers, verbose);
obj = neuralNetwork.convertDLToolboxNetwork(layers, verbose);

end


% Auxiliary functions -----------------------------------------------------

function dlt_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork)
    % read ONNX network via DLT

    % build name-value pairs
    NVpairs = {};

    % input data format
    if ~isempty(inputDataFormats)
        NVpairs = [NVpairs, {'InputDataFormats', inputDataFormats}];
    end

    % output data format
    if ~isempty(outputDataFormats)
        NVpairs = [NVpairs, {'OutputDataFormats', outputDataFormats}];
    end

    % custom layers generated from DLT will be stored in this folder
    % https://de.mathworks.com/help/deeplearning/ref/importnetworkfromonnx.html#mw_ccdf29c9-84cf-4175-a8ce-8e6ab1c89d4c
    customLayerName = 'DLT_CustomLayers';
    if isMATLABReleaseOlderThan('R2024a')
        % legacy
        NVpairs = [NVpairs, {'PackageName',customLayerName}];
    else
        NVpairs = [NVpairs, {'NameSpace',customLayerName}];
    end

    % load network
    if isMATLABReleaseOlderThan('R2023b')
        % legacy
        dlt_net = importONNXNetwork(file_path, NVpairs{:}, 'TargetNetwork', targetNetwork);
    else
        dlt_net = importNetworkFromONNX(file_path, NVpairs{:});
    end
end

function aux_removeIndentCodeLines(ME)

    % remove 'indentcode' function call 

    files = {'nnet.internal.cnn.onnx.fcn.ModelTranslation', 'nnet.internal.cnn.onnx.CustomLayerManager'};

    for i=1:length(files)
        % error happens in this file
        internalPath = which(files{i});
    
        % read text and comment failing line
        filetext = fileread(internalPath);
        filetext = strrep(filetext,"indentcode(","(");
    
        % try to read file with write permission
        fid  = fopen(internalPath,'w');
        if fid == -1
            % rethrowing error
            rethrow(ME)
        end
    
        % write new filetext
        fprintf(fid,'%s',filetext);
        fclose(fid);

    end

end

function layers = aux_groupCompositeLayers(layerslist, connections)
    % We use a digraph as a utility to convert the table of connections to
    % a nested cell array.

    % Build weights; we encode the order of inputs through weights in the
    % digraph.
    ts = regexp(connections.Destination, '/in(\d+)', 'tokens');
    weights = ones(length(ts), 1);
    % Logical index for non-empty tokens
    nonEmptyIdx = ~cellfun(@isempty, ts);
    weights(nonEmptyIdx) = cellfun(@(x) str2double(x{1}{1}), ts(nonEmptyIdx));

    % Strip the names. 
    srcs = cellfun(@(src) cell2mat(src), ...
        regexp(connections.Source,'^[^/]+','match'), ...
            'UniformOutput',false);
    dests = cellfun(@(dest) cell2mat(dest), ...
        regexp(connections.Destination,'^[^/]+','match'), ...
            'UniformOutput',false);

    % Convert to digraph.
    nodesTable = table({layerslist.Name}','VariableNames', {'Name'});
    G = digraph(srcs,dests,weights,nodesTable);
        
    % Sort the nodes topologically to find the final node.
    N = toposort(G);
    % Find the final node.
    nK = N(end);

    % % Visualization for debugging.
    % figure; 
    % plot(G,'EdgeLabel', G.Edges.Weight, ...
    %     'NodeLabel', cellfun(@(name) findnode(G,name),G.Nodes.Name));

    % Use a reverse depth-first traversal to construct a nested cell array.
    [layers,~] = aux_buildNestedCellArray(num2cell(layerslist), ...
        G,nK,[]);
end

function [layers,G,ni,currEdgeIdx] = aux_buildNestedCellArray(dltLayers,G, ...
    nK,currEdgeIdx)
    % We use a reverse depth-first traversal to construct a nested cell
    % array. We have to reverse the traversal to avoid not creating nested 
    % cells. 

    % Initialize the cell array.
    layers = {};
    % Initialize with last node.
    ni = nK;
    % Count the number of loop iteration to prevent an infinite loop.
    iter = 1;
    % Traverse the graph and build the nested cell array.
    while iter < height(G.Nodes)

        if ~isempty(currEdgeIdx)
            % Mark the current edge as visited by setting the weight to 0.
            G.Edges.Weight(currEdgeIdx) = 0;
            
            % % Visualization for debugging.
            % plot(G,'EdgeLabel', G.Edges.Weight, ...
            %     'NodeLabel', cellfun(@(name) findnode(G,name),G.Nodes.Name));
        end

        % Obtain the successors.
        outEdgeIdx = outedges(G,ni);

        if length(outEdgeIdx) > 1 
            % There are multiple computation paths starting at this node; 
            % we reached a forking node.
            if ~isempty(setdiff(outEdgeIdx,currEdgeIdx)) % any(G.Edges.Weight(outEdgeIdx) > 0)
                % There are still outgoing computation paths that have 
                % not been visited. We break out of the current 
                % computation path.
                break;
            end
        end

        % Prepend the current layer.
        layers = [dltLayers(ni) layers];

        % Obtain all incoming edges.
        inEdgeIdx = inedges(G,ni);

        if isempty(inEdgeIdx)
            % There are no incoming edges; we have reached the input node.
            break;
        end

        if length(inEdgeIdx) > 1
            % There are multiple computation paths merging in the current
            % node. We have to construct the layer cell arrays for each
            % computation path.
            layersPreds = {};

            % Find the weight of the edges, which encode the order 
            % (required for concatenation).
            ws = G.Edges.Weight(inEdgeIdx)';
            % Sort the predecessors based on the edge weight.
            [~,predIdx] = sort(ws,'ascend');

            % Store edges of each computation path.
            currEdgeIdx = [];

            % Create a new computation path for each predecessor.
            for j=predIdx
                % Mark the edge as visited by setting the weight to 0.
                G.Edges.Weight(inEdgeIdx(j)) = 0;
                % Obtain the start of the j-th predecessor edge.
                predj = findnode(G,G.Edges.EndNodes{inEdgeIdx(j),1});
                % Check if the computation path is just a residual
                % connection (has multiple outedges).
                isResidualConn = length(outedges(G,predj)) > 1;
                if isResidualConn
                    % Prepend the empty cell for the j-th computation path.
                    layersPreds = [layersPreds; {{}}];
                    % Append the residual edge.
                    currEdgeIdx = [currEdgeIdx inEdgeIdx(j)];
                    % Update the current node.
                    nij = predj;
                else
                    % Compute the layers of the j-th computation path.
                    [layersj,G,nij,currEdgeIdxj] = ...
                        aux_buildNestedCellArray(dltLayers,G,predj,[]);
                    % Prepend the layers of the j-th computation path.
                    layersPreds = [layersPreds; {layersj}];
                    % Append final node.
                    currEdgeIdx = [currEdgeIdx currEdgeIdxj];
                end
            end
            % Traversed all computation paths; prepend the layers.
            layers = [{layersPreds} layers];
            % Update the current node.
            ni = nij;
        else
            % There is only a single incoming edge; move to the only 
            % predecessor.
            currEdgeIdx = inEdgeIdx(1);
            ni = findnode(G,G.Edges.EndNodes{inEdgeIdx(1),1});
        end

        % Increment iteration counter.
        iter = iter + 1;
    end
end

% ------------------------------ END OF CODE ------------------------------
