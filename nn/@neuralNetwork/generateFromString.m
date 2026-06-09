function nn = generateFromString(layersstring,inSz,nK,varargin)
% generateFromString - creates a random layer-based network
%
% Syntax:
%    % Specify the layers.
%    layersstring = [
%        'CONV 5 4*4+2' % conv. layer with 5 filters of size 4x4 and stride 2
%        'CONV 10 4*4+1'
%        'FC 100' % fully connected layer with 100 neurons
%    ];
%    batchnorm = true; % Add batchnorm before every activation function.
%    actfun = 'relu'; % Specify type of activation function.
%    % Generate the neural network.
%    nn = neuralNetwork.generateFromString(layertypes,[28 28 1]);
%
% Inputs:
%   - layersstring: string where each line specifies a layer type, i.e.,
%    * 'FC n' for fully connected layers with n neurons
%    * 'CONV k w*h+s' for conv. layer with k filters of size w*h and
%       stride s
%   - inSz: array specifying the size of an input image,
%    e.g., [28 28 1] for MNIST.
%   - nK: number of output dimensions.
%   - actfun: type of activation function that is added after each layer,
%       i.e., {'relu','tanh',...}.
%   - batchnorm: add a batch norm before each activation function.
%   - verbose: verbose output.
%
% Outputs:
%    nn - generated neural network
%    layers - the layers of the neural network
%    nn - generated neural network
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork, nnActivationLayer/instantiateFromString

% Authors:       Lukas Koller
% Written:       17-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Set default arguments.
[actfun,batchnorm,verbose] = setDefaultValues({'relu',false,false},varargin);

% Split the given string to a cell array, each entry corresponds to a line.
layertypes = strtrim(cellstr(splitlines(string(layersstring))));
% Remove empty lines.
layertypes = layertypes(~cellfun(@isempty,layertypes));

% Specify regex patterns to match the parameters of the layer types.
fcExpr = '(?<=FC\s)\d+';
convExpr = 'CONV\s+(?<k>\d+)\s+(?<w>\d+)\*(?<h>\d+)\+(?<s>\d+)';

% Initialize a cell array to store the layers.
layers = {};
% Initialize the size of the feature map.
n = inSz;

% Instantiate the different layer types.
for i = 1:length(layertypes)
    % Obtain the current layer type.
    typei = layertypes{i};

    % Check the type of the current layer.
    if startsWith(typei, 'FC')
        % Extract the number of neurons.
        ni = str2double(regexp(typei,fcExpr,'match','once'));
        % Initialize zero weights and bias.
        Wi = zeros([ni prod(n)]);
        bi = zeros([ni 1]);
        % Create and append a linear layer.
        layers{end+1} = nnLinearLayer(Wi,bi);
        % Update the feature map size.
        n = [ni 1];
    elseif startsWith(typei, 'CONV')
        % Extract the parameters for the convolutional layer, i.e.,
        % number of filter k, kernel size w*h, and stride s.
        psi = regexp(typei,convExpr,'names');
        % Convert the parameters.
        ki = str2double(psi.k);
        wi = str2double(psi.w);
        hi = str2double(psi.h);
        si = str2double(psi.s);
        % Initialize zero filters and bias.
        Wi = zeros([wi hi n(end) ki]);
        bi = zeros([ki 1]);
        % Create and append a convolutional layer.
        layers{end+1} = nnConv2DLayer(Wi,bi,[0 0 0 0],[si si]);
        % Update the feature map size.
        n = layers{end}.getOutputSize(n);
    end

    % Append a batch normalization layer.
    if batchnorm
        layers{end+1} = nnBatchNormLayer();
    end

    % Append and append an activation layer.
    layers{end+1} = nnActivationLayer.instantiateFromString(actfun);
end

% Initialize zero weights and bias for the output layer.
WK = zeros([nK prod(n)]);
bK = zeros([nK 1]);
% Add the final output layer.
layers{end+1} = nnLinearLayer(WK,bK);

% Instantiate the neural network.
nn = neuralNetwork(layers);
% Set the input size of the neural network.
nn.setInputSize(inSz,verbose);
% Initialize the weights and bias.
nn.initWeights('glorot');

% ------------------------------ END OF CODE ------------------------------
