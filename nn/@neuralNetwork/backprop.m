function grad_in = backprop(nn, grad_out, options, varargin)
% backprop - compute the backpropagation for the previous input
%
% Syntax:
%    grad_in = nn.backprop(grad_out, options, idxLayer)
%
% Inputs:
%    nn - neuralNetwork
%    grad_out - gradient of the output of the neural network
%    options - training parameters
%    idxLayer - indices of layers that should be evaluated
%
% Outputs:
%    grad_in - gradient w.r.t the input
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate, nnOptimizer

% Authors:       Tobias Ladner, Lukas Koller
% Written:       01-March-2023
% Last update:   03-May-2023 (LK, added backprop for polyZonotope)
%                25-May-2023 (LK, added options as function parameter)
%                31-July-2023 (LK, return update gradients)
%                04-August-2023 (LK, added layer indices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,4)
if nargin == 2
    options = struct;
end

inputArgsCheck({ ...
    {nn, 'att', 'neuralNetwork'}; ...
    {grad_out, 'att', {'numeric','interval','gpuArray'}}
    {options, 'att', 'struct'}; ... 
})

% validate parameters
[idxLayer] = setDefaultValues({1:length(nn.layers)}, varargin);

% execute -----------------------------------------------------------------

if isnumeric(grad_out)
    % numeric
    for i = flip(idxLayer)
        layer_i = nn.layers{i};
        % Retrieve stored input
        if options.nn.train.backprop
            input = layer_i.backprop.store.input;
        end
        grad_out = layer_i.backpropNumeric(input, grad_out, options);
    end
    grad_in = grad_out;
else
    throw(CORAerror('CORA:notSupported',...
        ['Set representation ' class(grad_out) ' is not supported.']));
end

end

% ------------------------------ END OF CODE ------------------------------
