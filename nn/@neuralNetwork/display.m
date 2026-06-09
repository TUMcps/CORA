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
% Last update:   17-January-2023
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
aux_printLayers(obj.layers, 0);
fprintf("\n")

end


% Auxiliary functions -----------------------------------------------------

function aux_printLayers(layers, indentLevel)
    for i = 1:length(layers)
        layer_i = layers{i};
        
        % create indentation string
        indent = repmat('    ', 1, indentLevel);

        fprintf("%s (%d)\t %s\n", indent, i, layer_i.getLayerInfo())
        
        % check for composite layer
        if isa(layer_i, 'nnCompositeLayer')
            for j = 1:length(layer_i.layers)
                fprintf("%s   Path %d:\n", indent, j);
                aux_printLayers(layer_i.layers{j}, indentLevel + 1);
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
