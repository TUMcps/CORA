function nn_normal = getNormalForm(obj)
% getNormalForm - transforms the neural network into an identical network 
%    in normal form by combining subsequent linear layers s.t. there are 
%    only alternating linear and nonlinear layers.
%    Note: only works on feed-forward networks
%
% Syntax:
%    nn_normal = neuralNetwork.getNormalForm(obj)
%
% Inputs:
%    obj - neuralNetwork
%
% Outputs:
%    nn_normal - neural network in normal form
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       14-December-2022
% Last update:   17-January-2023 (TL, Reshape)
% Last revision: 01-August-2023

% ------------------------------ BEGIN CODE -------------------------------

layers = {};

% init properties
W = 1; 
b = 0;

for i=1:length(obj.layers)
    layer = obj.layers{i};

    if isa(layer, 'nnConv2DLayer')
        layer = layer.convert2nnLinearLayer();
    end

    if isa(layer, 'nnLinearLayer')
        W = layer.W * W;
        b = layer.W * b;
        b = sum(b, 2); % required to make b a column vector
        b = b + layer.b;
    elseif isa(layer, 'nnElementwiseAffineLayer')
        W = diag(layer.scale) * W;
        b = layer.scale .* b + layer.offset;
    elseif isa(layer, 'nnIdentityLayer')
        % W = W; b = b;
    elseif isa(layer, 'nnReshapeLayer')
        W_reshape = eye(prod(layer.inputSize));
        W_reshape = layer.evaluateNumeric(W_reshape, struct);
        W = W_reshape * W;

    elseif isa(layer, 'nnActivationLayer')
        % create new linear layer
        layers{end+1} = aux_getLinearLayer(W, b);

        % reset properties
        W = 1; b = 0;

        % activation
        layers{end+1} = layer.copy();
    else
        throw(CORAerror('CORA:notSupported', ...
            sprintf("Unable to bring a layer of type '%s' in normal form.", class(layer))))
    end
end

layer = aux_getLinearLayer(W, b);
if ~isa(layer, 'nnIdentityLayer')
    layers{end+1} = layer;
end

% construct neural network in normal form
nn_normal = neuralNetwork(layers);

end


% Auxiliary functions -----------------------------------------------------

function layer = aux_getLinearLayer(W, b)
    if length(W) > 1 || length(b) > 1
        layer = nnLinearLayer(W, b);
    elseif W ~= 1 || b ~= 0
        layer = nnElementwiseAffineLayer(W, b);
    else
        layer = nnIdentityLayer();
    end
end

% ------------------------------ END OF CODE ------------------------------
