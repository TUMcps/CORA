function numNeurons = getNumNeurons(obj)
% getNumNeurons - returns the number of input neurons per layer and the
%    number of output neurons
%
% Syntax:
%    pattern = getNumNeurons(obj)
%
% Inputs:
%    obj - object of class neuralNetwork
%
% Outputs:
%    numNeurons - array of number of neurons per layer, or nan if unknown
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-August-2022
% Last update:   ---
% Last revision: 02-October-2023

% ------------------------------ BEGIN CODE -------------------------------

numNeurons = zeros(1, size(obj.layers, 1)+1);

% iterate over all layers
for i = 1:length(obj.layers)    
    if isempty(obj.layers{i}.inputSize)
        % try to set input size for all layers
        try
            obj.setInputSize();
        end
    end

    if isempty(obj.layers{i}.inputSize)
        % if still empty, use nan
        numNeurons(i) = nan;
    else
        numNeurons(i) = prod(obj.layers{i}.inputSize);
    end
end

% output neurons
if isempty(obj.neurons_out)
    numNeurons(end) = nan;
else
    numNeurons(end) = obj.neurons_out;
end

end

% ------------------------------ END OF CODE ------------------------------
