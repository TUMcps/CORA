function propagateBounds(obj, i, varargin)
% propagateBounds - propagates the bounds of layer i to the next activation
%   layer and updates the bounds accordingly
%
% Syntax:
%    neuralNetwork.propagateBounds(obj, i)
%    neuralNetwork.propagateBounds(obj, i, evParams)
%
% Inputs:
%    obj - neuralNetwork
%    i - starting layer
%    evParams - struct, evaluation parameters
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluate

% Authors:       Tobias Ladner
% Written:       15-December-2022
% Last update:   ---
%                ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs', 3))
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs', 4))
end
[evParams] = setDefaultValues({struct}, varargin);
inputArgsCheck({ ...
    {obj, 'att', 'neuralNetwork'}; ...
    {i, 'att', 'numeric', 'scalar'}; ...
    {evParams, 'att', 'struct'}; ...
})

layer_i = obj.layers{i};
if evParams.propagate_bounds && evParams.reuse_bounds ...
    && isa(layer_i, 'nnActivationLayer')
    l = layer_i.l;
    u = layer_i.u;

    if ~isempty([l; u]) && ~any(isnan([l; u]))

        bounds = interval(l, u);
        bounds = layer_i.evaluateInterval(bounds, evParams);

        % propagate bounds
        for j = i + 1:length(obj.layers)
            layer_j = obj.layers{j};

            if isa(layer_j, 'nnActivationLayer')
                if isempty(layer_j.l) || any(isnan(layer_j.l)) || ...
                        isempty(layer_j.u) || any(isnan(layer_j.u))
                    % set bounds
                    layer_j.l = bounds.inf;
                    layer_j.u = bounds.sup;

                else
                    % update bounds
                    Ilayer = interval(layer_j.l, layer_j.u);
                    bounds = bounds & Ilayer;
                    layer_j.l = bounds.inf;
                    layer_j.u = bounds.sup;
                end
                return
            end

            bounds = layer_j.evaluateInterval(bounds, evParams);
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
