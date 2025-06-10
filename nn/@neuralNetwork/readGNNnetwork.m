function obj = readGNNnetwork(file_path, varargin)
% readGNNnetwork - reads and converts a gnn network into
% a cora neuralNetwork for verification
%
% Syntax:
%    obj = neuralNetwork.readGNNnetwork(file_path,verbose)
%
% Inputs:
%    file_path - path to file (list of layer objects)
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
%                29-January-2024 (TL, new model structure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,2);
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
layer_number = length(model.layers);
layers = {};

% iterate through the fields of the file to check and append layers together
for i = 1:layer_number

    % read data
    layer_i = model.layers(i);
    type_i = layer_i.type;
    W_i = layer_i.W;
    b_i = layer_i.b;
    act_i = layer_i.act;

    % linear layer
    switch(type_i)
        case "gcn"
            layers{end+1} = nnGCNLayer();
            layers{end+1} = nnGNNLinearLayer(W_i, b_i);

        case "lin"
            layers{end+1} = nnGNNLinearLayer(W_i, b_i);

        case 'global_add_pool'
            layers{end+1} = nnGNNGlobalPoolingLayer('add');

        case 'global_mean_pool'
            layers{end+1} = nnGNNGlobalPoolingLayer('mean');

        otherwise
            throw(CORAerror('CORA:wrongFieldValue', type_i, {'gcn', 'lin','global_add_pool','global_mean_pool'}))
    end

    % following by activation
    if ~isempty(act_i)
        layers{end+1} = nnActivationLayer.instantiateFromString(act_i);
    end
end

obj = neuralNetwork(layers);

% print the layers if verbose is true
if verbose
    display(obj)
end

end

% ------------------------------ END OF CODE ------------------------------
