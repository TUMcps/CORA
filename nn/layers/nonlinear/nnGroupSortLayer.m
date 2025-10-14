classdef nnGroupSortLayer < nnActivationLayer
% nnGroupSortLayer - class for GroupSort layers. For grSz = 2 it is
%   equivalent to a MinMaxLayer. For grSz=inputSize, equivalent to 
%   a FullSortLayer.
%
% Syntax:
%    obj = nnGroupSortLayer(grSz,name)
%
% Inputs:
%    grSz - size of the sorted groups/chunks.
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Anil, C. et al. Sorting out Lipschitz function approximation. (ICML). 2019.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Daniel Safyan, Lukas Koller
% Written:       18-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    groupSize = 2       % size of the sorted chunks 
end

methods
    function obj = nnGroupSortLayer(varargin)
        % Parse input arguments.
        [groupSize, name] = setDefaultValues({2, []}, varargin);

        % Call super class constructor.
        obj@nnActivationLayer(name)
        obj.groupSize = groupSize;
    end

    function df_i = getDf(obj, i)
        % The activation function is not elementwise; thus, not well
        % defined.
        df_i = 1;
    end

    function [df_l,df_u] = getDerBounds(obj, l, u)
        df_l = 1;
        df_u = 1;
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    function r = evaluateNumeric(obj, input, options)
        % Obtain the number of input dimensions and batch size.
        [n,bSz] = size(input);
        % Compute the number of elements which do not fit in the last
        % group.
        numGroups = floor(n/obj.groupSize);
        numRes = mod(n,obj.groupSize);
        numPad = obj.groupSize - numRes;
        % Obtain padded number of dimensions.
        n_ = n + numPad;
        % Pad the input to fill up the last group.
        input_ = [input; Inf([numPad bSz])];
        % Compute the number of groups.
        numGroups_ = n_/obj.groupSize;
        % Reshape the input s.t. each group is located in one dimension.
        input_ = reshape(input_,[obj.groupSize numGroups_ bSz]);
        % Sort along the group dimension.
        [r_,idx_] = sort(input_,1,'ascend');
        % Reshape the sorted output to correct shape.
        r = reshape(r_,[n_ bSz]);
        % Trim the padded dimensions.
        r = r(1:n,:);
        
        if options.nn.train.backprop
            % Compute permutation indices.
            idx_ = sub2ind([obj.groupSize numGroups_ bSz],idx_, ...
                repelem(1:numGroups_,obj.groupSize,1,bSz), ...
                permute(repelem(1:bSz,obj.groupSize,1,numGroups_),[1 3 2]));
            % We convert the sort indices to indices into input. Therefore,
            % we construct a dummy input which contains indices instead of
            % values.
            inputIdx = reshape( ...
                [reshape(1:numel(input),size(input)); Inf([numPad bSz])], ...
                [obj.groupSize numGroups_ bSz]);
            % Reshape into the correct dimensions.
            idx = reshape(inputIdx(idx_),[n_ bSz]);
            % Trim the padded dimensions.
            idx = idx(1:n,:);
            % Store the permutation indices for backpropagation.
            obj.backprop.store.permIdx = idx;
        end
    end


% backprop ----------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options, updateWeights)
        % Obtain the stored permutation.
        permIdx = obj.backprop.store.permIdx;
        % Apply the permutation to the gradient.
        grad_in = grad_out(permIdx);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
