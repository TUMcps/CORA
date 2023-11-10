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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% validate parameters
if nargin < 1
    throw(CORAerror("CORA:notEnoughInputArgs", 1))
elseif nargin > 5
    throw(CORAerror("CORA:tooManyInputArgs", 5))
end
[verbose, inputDataFormats, outputDataFormats, targetNetwork] = ...
    setDefaultValues({false, 'BC', 'BC', 'dagnetwork'}, varargin);

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

% convert DLT network to CORA network
obj = neuralNetwork.convertDLToolboxNetwork(dltoolbox_net.Layers, verbose);


end


% Auxiliary functions -----------------------------------------------------

function dltoolbox_net = aux_readONNXviaDLT(file_path,inputDataFormats,outputDataFormats,targetNetwork)
    % read ONNX network via DLT

    % custom layers generated from DLT will be stored in this folder
    % https://de.mathworks.com/help/deeplearning/ref/importnetworkfromonnx.html#mw_ccdf29c9-84cf-4175-a8ce-8e6ab1c89d4c
    customLayerName = 'DLT_CustomLayers';

    if isempty(which('importNetworkFromONNX'))
        % legacy
        dltoolbox_net = importONNXNetwork(file_path, ...
            'InputDataFormats', inputDataFormats, ...
            'OutputDataFormats', outputDataFormats, ...
            'PackageName', customLayerName, ...
            'TargetNetwork', targetNetwork);
    else
        % MATLAB >=R2023b
        dltoolbox_net = importNetworkFromONNX(file_path, ...
            'InputDataFormats', inputDataFormats, ...
            'OutputDataFormats', outputDataFormats, ...
            'PackageName', customLayerName);
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

% ------------------------------ END OF CODE ------------------------------
