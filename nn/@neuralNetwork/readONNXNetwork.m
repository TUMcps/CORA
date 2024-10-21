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

% Authors:       Tobias Ladner
% Written:       30-March-2022
% Last update:   07-June-2022 (specify in- & outputDataFormats)
%                30-November-2022 (removed neuralNetworkOld)
%                13-February-2023 (simplified function)
%                21-October-2024 (clean up DLT function call)
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
    dltoolbox_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork);

catch ME
    if strcmp(ME.identifier, 'MATLAB:javachk:thisFeatureNotAvailable') && ...
            contains(ME.message,'Swing is not currently available.')
        % matlab tries to indent the code of the generated files for 
        % custom layers, for which (somehow?) a gui is required.
        % As e.g. docker runs don't have a gui, we just try to remove the
        % 'indentcode' function call here ...
        aux_removeIndentCodeLines(ME);

        % re-read network
        dltoolbox_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork);

    else
        rethrow(ME)
    end
end

if containsCompositeLayers
    % Combine multiple layers into blocks to realize residual connections and
    % parallel computing paths.
    layers = aux_groupCompositeLayers(dltoolbox_net.Layers,dltoolbox_net.Connections);
else
    layers = num2cell(dltoolbox_net.Layers);
end

% convert DLT network to CORA network
% obj = neuralNetwork.convertDLToolboxNetwork(dltoolbox_net.Layers, verbose);
obj = neuralNetwork.convertDLToolboxNetwork(layers, verbose);


end


% Auxiliary functions -----------------------------------------------------

function dltoolbox_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork)
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
        dltoolbox_net = importONNXNetwork(file_path, NVpairs{:}, 'TargetNetwork', targetNetwork);
    else
        dltoolbox_net = importNetworkFromONNX(file_path, NVpairs{:});
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
    % Find initial layer.
    layer0 = aux_findLayerByName(layerslist,connections(1,:).Source{1});
    layers = {layer0};
    for i=1:height(connections)
        % Find source and destination layer.
        % layerSrc = aux_findLayerByName(layerslist,connections(i,:).Source{1});
        layerDest = aux_findLayerByName(layerslist,connections(i,:).Destination{1});
        % Check if the next layer aggreates multiple paths.
        isAggrLayer = ~strcmp(layerDest.Name,connections(i,:).Destination{1});
        % Find source layer within the current list of layers.
        for j=length(layers):-1:1 % Speed up computations by searching from the back.
            % Obtain paths from the j-th step.
            layerj = layers{j};
            if iscell(layerj)
                % Iterate over the current paths to try and find the source.
                for k=1:length(layerj)
                    % Extract layer from the j-th step and k-th path.
                    layerjk = layerj{k};
                    for l=length(layerjk):-1:1
                        if strcmp(layerjk(l).Name,connections(i,:).Source{1})
                            % Append in computation path.
                            if isAggrLayer
                                % Add at the end.
                                layers{j+1} = layerDest;
                            else
                                layers{j}{k} = [layerjk layerDest];
                            end
                            break;
                        end
                    end
                end
            else
                % Compare the names.
                if strcmp(layerj.Name,connections(i,:).Source{1})
                    % Found source layer; append the destination.
                    if j+1 > length(layers)
                        % Add at the end.
                        layers{j+1} = layerDest;
                    else
                        if isAggrLayer
                            % There is a residual connection.
                            layers{j+1} = {layers{j+1}; []};
                            layers{j+2} = layerDest;
                        else
                            % There is a new computation path.
                            layers{j+1} = {layers{j+1}; layerDest};
                        end
                    end
                    break;
                elseif startsWith(layerj.Name,connections(i,:).Source{1})
                    % Combine computation paths.
                    % TODO
                    break;
                end
            end
            
            % if iscell(layerj)
            %     % For now we do not supported nested residual connections.
            %     continue;
            % end
        end
    end    
end

function layer = aux_findLayerByName(layers, destName)
    idx = regexp(destName,'/*','once');
    if ~isempty(idx)
        destName(idx:end) = [];
    end
    layer = [];
    for i=1:length(layers)
        % if startsWith(destName,layers(i).Name)
        if strcmp(destName,layers(i).Name) ... || regexp(destName,[layers(i).Name '/*'])
            layer = layers(i);
            break;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
