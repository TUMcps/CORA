function display(obj)
% display - displays the properties of a neuralNetwork object
%
% Syntax:
%    display(obj)
%
% Inputs:
%    obj - neuralNetwork object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       23-November-2022
% Last update:   17-January-2023 (TL, better layer output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% disp input if necessary
dispInput(inputname(1))

fprintf("Neural network: '%s'\n", class(obj))

% in/out neurons
fprintf("Nr. of input neurons: %d\n", obj.neurons_in)
fprintf("Nr. of output neurons: %d\n", obj.neurons_out)
fprintf(newline);

fprintf("layers: (%d layers)\n", length(obj.layers))
for i = 1:length(obj.layers)
    layer_i = obj.layers{i};
    fprintf(" (%d)\t %s\n", i, layer_i.getLayerInfo())
end
fprintf("\n")

% ------------------------------ END OF CODE ------------------------------
