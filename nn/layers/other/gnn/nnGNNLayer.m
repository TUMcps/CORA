classdef (Abstract) nnGNNLayer < nnLayer
% nnGNNLayer - abstract class for nn layers
%
% Syntax:
%    layer = nnGNNLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Ladner, T., et al. (2025). Formal Verification of Graph 
%        Convolutional Networks with Uncertain Node Features 
%        and Uncertain Graph Structure. TMLR.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_gnn_uncertain_message_passing

% Authors:       Gerild Pjetri, Tianze Huang, Tobias Ladner
% Written:       13-December-2022
% Last update:   23-February-2023 (TL, integrate in nnLayer)
% Last revision: 15-January-2023

% ------------------------------ BEGIN CODE -------------------------------

methods
    function obj = nnGNNLayer(varargin)
        obj@nnLayer(varargin{:})
        obj.inputSize = [];
    end
end

methods (Abstract)
    [nin, nout] = getNumNeurons(obj, graph)
    outputSize = getOutputSize(obj, inputSize, graph)
end

end

% ------------------------------ END OF CODE ------------------------------
