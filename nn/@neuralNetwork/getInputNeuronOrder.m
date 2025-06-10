function neuronOrder = getInputNeuronOrder(nn,method,x,varargin)
% getInputNeuronOrder - ranks input neurons based on the specified method
%
% Syntax:
%    neuronOrder = getInputNeuronOrder(nn,method,x)
%
% Inputs:
%    nn - object of class neuralNetwork
%    method - str, method to determine neuron order: 
%        'in-order', 'sensitivity', 'snake'
%    x - numeric vector
%    inputSize - numeric 
%
% Outputs:
%    -
% 
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       04-July-2024
% Last update:   15-November-2024 (added snake)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(3,4)
[inputSize] = setDefaultValues({[]},varargin);
inputArgsCheck({ ...
    {nn,'att','neuralNetwork'}; ...
    {method,'str',{'in-order','sensitivity','snake'}}; ...
    {x,'att','numeric'}; ...
    {inputSize,'att','numeric'}; ...
})

% determine order based on method
switch method

    case 'in-order'
        % just iterate through all pixels
        if isempty(inputSize)
            neuronOrder = 1:numel(x);
        else
            neuronOrder = 1:(inputSize(1)*inputSize(2));
        end

    case 'snake'
        % spiral inwards
        featOrderNaive = 1:(inputSize(1)*inputSize(2));
        A = reshape(featOrderNaive,inputSize(1:2));

        neuronOrder = [];
        while ~isempty(A)
           neuronOrder = [neuronOrder,A(1,:)];
           A = fliplr(A(2:end,:))';
        end        

    case 'sensitivity' 
        % identify least sensitive input neurons
        S = nn.calcSensitivity(x);

        % mean across output neurons
        S = mean(abs(S),1); 
        if ~isempty(inputSize)
            % mean across channels
            S = reshape(S,inputSize);
            S = mean(S,3);
            S = reshape(S,[],1);
        end
        [~,neuronOrder] = sort(S);

    otherwise
        % should not happen as it was checked before
        throw(CORAerror('CORA:specialError',sprintf('Unknown method: %s', featOrderMethod)))

end    

end

% ------------------------------ END OF CODE ------------------------------
