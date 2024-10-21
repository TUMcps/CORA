function propagateBounds(obj, i, varargin)
% propagateBounds - propagates the bounds of layer i to the next activation
%   layer and updates the bounds accordingly
%
% Syntax:
%    neuralNetwork.propagateBounds(obj, i)
%    neuralNetwork.propagateBounds(obj, i, options)
%
% Inputs:
%    obj - neuralNetwork
%    i - starting layer
%    options - struct, evaluation parameters (stored in options.nn)
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
narginchk(3,4);
[options] = setDefaultValues({struct}, varargin);
inputArgsCheck({ ...
    {obj, 'att', 'neuralNetwork'}; ...
    {i, 'att', 'numeric', 'scalar'}; ...
    {options, 'att', 'struct'}; ...
})

layer_i = obj.layers{i};
if options.nn.propagate_bounds && options.nn.reuse_bounds ...
    && isa(layer_i, 'nnActivationLayer')
    l = layer_i.l;
    u = layer_i.u;

    if ~isempty([l; u]) && ~any(isnan([l; u]))

        bounds = interval(l, u);
        bounds = layer_i.evaluateInterval(bounds, options);

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

            bounds = layer_j.evaluateInterval(bounds, options);
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
