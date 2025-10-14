function initWeights(nn,varargin)
% initWeights - initialize the weight of the neural network
%
% Syntax:
%    nn.initWeights('glorot')
%
% Inputs:
%    method - {'glorot' (default), 'shi'}
%    seed - 'default' or nonnegative integer
%    idxLayer - indices of layers that should be initialized
%
% Outputs:
%    - 
%
% References:
%    [1] Glorot, X., Bengio, Y. "Understanding the difficulty of training 
%       deep feedforward neural networks". PLMR, 2010.
%    [2] Z. Shi, Y. Wang, H. Zhang, J. Yi, and C. Hsieh, "Fast certified 
%       robust training with short warmup". NeurIPS, 2021.
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/train

% Authors:       Lukas Koller
% Written:       12-February-2024
% Last update:   18-August-2025 (enumerate layers)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Enumerate all layers.
[layers,~] = nn.enumerateLayers();

% validate parameters
[method,seed,idxLayer] = setDefaultValues({'glorot','default',...
    1:length(layers)}, varargin);

rng(seed); % Set seed for reproducability.

for i=idxLayer
    layeri = layers{i};
    if isempty(layeri.getLearnableParamNames())
        % There are no learnable parameters to initialize.
        continue;
    end

    if isa(layeri,'nnLinearLayer') ...
        || isa(layeri,'nnConv2DLayer') ...
        || isa(layeri,'nnLipConstrLinearLayer')
        [nin, ~] = layeri.getNumNeurons();
        if strcmp(method,'glorot')
            % uniform between [-a,a] where a = 1/sqrt(nin)
            a = 1/sqrt(nin);
            layeri.W = unifrnd(-a,a,size(layeri.W));
            % init bias with 0
            layeri.b = zeros(size(layeri.b));
        elseif strcmp(method,'shi')
            % normal distributed with mu = 0, sigma = sqrt(2*pi)/nin
            sigma = sqrt(2*pi)/nin;
            layeri.W = normrnd(0,sigma,size(layeri.W));
            % init bias with 0
            layeri.b = zeros(size(layeri.b));
        end
    elseif isa(layeri,'nnBatchNormLayer')
        layeri.scale = 1;
        layeri.offset = 0;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
