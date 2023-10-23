function obj = readGNNnetwork(file_path, varargin)
% readGNNnetwork - reads and converts a gnn network into
% a cora neuralNetwork for verification
%
% Syntax:
%    obj = readGNNnetwork(file_path,verbose)
%
% Inputs:
%    file_path - path to file(information of # of layers, bias, weight...'features')
%    verbose - bool if information should be displayed
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/readONNXNetwork

% Authors:       Gerild Pjetri, Tianze Huang, Tobias Ladner
% Written:       02-December-2022
% Last update:   15-January-2023
%                23-February-2023 (TL, clean up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs', 1))
elseif nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs', 2))
end
verbose = setDefaultValues({false}, varargin);
inputArgsCheck({ ...
    {file_path, 'att', {'char', 'string'}}; ...
    {verbose, 'att', 'logical'}; ...
})
if ~isfile(file_path)
    throw(CORAerror('CORA:fileNotFound', file_path));
end

% prepare and setup the model file containing weights and biases
fid = fopen(file_path);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
model = jsondecode(str);

% layers: cell array containing different types of layers

% calculating the number of layers
fn = fieldnames(model);
l = length(fn);
layer_number = (l - 4) / 4;
layers = {};

% iterate through the fields of the file to check and append layers together
for i = 0:layer_number
    type_name = aux_getJSONLayerName(i, 'type');
    type_i = model.(type_name);
    % linear layer
    if type_i == "gcn"
        W = double(model.(aux_getJSONLayerName(i, 'weight')));
        b = double(model.(aux_getJSONLayerName(i, 'bias')));
        layers{end+1} = nnGCNLayer(W, b);

    elseif type_i == "lin"
        W = double(model.(aux_getJSONLayerName(i, 'weight')));
        b = double(model.(aux_getJSONLayerName(i, 'bias')));
        layers{end+1} = nnGNNLinearLayer(W, b);

    else
        throw(CORAerror('CORA:wrongFieldValue', type_name, {'gcn', 'lin'}))
    end

    % following by activation
    act_name = aux_getJSONLayerName(i, 'act');
    if isfield(model, act_name)
        act_i = model.(act_name);
        layers{end+1} = nnActivationLayer.instantiateFromString(act_i);
    end
end

% append global pooling layer last
if isfield(model, 'global_pooling_type')
    layers{end+1} = nnGNNGlobalPoolingLayer(model.global_pooling_type);
end

obj = neuralNetwork(layers);

% print the layers if verbose is true
if verbose
    display(obj)
end

end


% Auxiliary functions -----------------------------------------------------

function str = aux_getJSONLayerName(n, type)
    str = ['l_', num2str(n), '_', type];
end

% ------------------------------ END OF CODE ------------------------------
