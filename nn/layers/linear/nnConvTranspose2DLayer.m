classdef nnConvTranspose2DLayer < nnConv2DLayer
% nnConvTranspose2DLayer - class for convolutional transpose 2D layers
%
% Inherits from nnConv2DLayer and swaps the forward/backward operations
% and overrides the output size calculation.
%
% Syntax:
%    obj = nnConvTranspose2DLayer(W, b, padding, stride, dilation, name)
%
% Inputs:
%    W - weight matrix (4-D single)
%        (kernel_height, kernel_width, out_channels, in_channels)
%    b - bias column vector (out_channels)
%    padding - zero padding [left top right bottom]
%    stride - step size per dimension
%    dilation - step size per dimension
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%
% Other m-files required: nnConv2DLayer
% Subfunctions: none
% MAT-files required: none

% See also: nnConv2DLayer

% Authors:       Benedikt Kellner
% Written:       01-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods
    function obj = nnConvTranspose2DLayer(varargin)
        if nargin >= 1 && isnumeric(varargin{1})
            % Permute W from (H,W,Out,In) to (H,W,In,Out) for parent
            W = varargin{1};
            varargin{1} = permute(W, [1 2 4 3]);
        end
        obj@nnConv2DLayer(varargin{:})
    end
end

methods (Access = protected)
    function [out_h, out_w, out_c] = aux_computeOutputSize(obj, varargin)
        [Filter, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({obj.W, obj.inputSize, obj.stride, ...
                obj.padding, obj.dilation}, varargin);
        
        in_h = inImgSize(1);
        in_w = inImgSize(2);
        
        [f_h, f_w] = obj.aux_getFilterSize(Filter, dilation);
        
        % padding [left,top,right,bottom]
        pad_l = padding(1);
        pad_t = padding(2);
        pad_r = padding(3);
        pad_b = padding(4);
        
        stride_h = stride(1);
        stride_w = stride(2);
        
        % Transposed convolution output size
        out_h = (in_h - 1)*stride_h + f_h - (pad_t + pad_b);
        out_w = (in_w - 1)*stride_w + f_w - (pad_l + pad_r);
        out_c = size(Filter, 4);
    end

    function [r,Wff] = conv2d(obj,input,options,varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ... 
                obj.stride, obj.padding, obj.dilation}, varargin);
        
        if options.nn.use_dlconv
            [~,batchSize] = size(input);
    
            % padding [left,top,right,bottom]
            pad_l = padding(1);
            pad_t = padding(2);
            pad_r = padding(3);
            pad_b = padding(4);
    
            inputImg = dlarray(reshape(input,[inImgSize batchSize]),'SSCB');
            if isempty(b)
                b = 0;
            end
            
            % Permute Filter: (H,W,In,Out) -> (H,W,Out,In) to match dltranspconv
            Filter = permute(Filter, [1 2 4 3]);

            rImg = dltranspconv(inputImg,Filter,b, ...
                Stride=stride,DilationFactor=dilation,...
                    Cropping=[pad_t pad_l; pad_b pad_r]);
            r = reshape(extractdata(rImg),[],batchSize);
    
            Wff = [];
        else
            throw(CORAerror('CORA:notImplemented', ...
                'Non-dlconv execution not supported yet.'));
        end
    end

    function r = transconv2d(obj,input,options,varargin)
        % Backward pass of transposed convolution, implemented as a
        % standard convolution (dlconv) with permuted weights.

        % Default input size is the layer's output size (the "large" side)
        defaultInSize = obj.getOutputSize(obj.inputSize);

        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, defaultInSize, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

        if options.nn.use_dlconv
            [~,batchSize] = size(input);

            inputImg = dlarray(reshape(input,[inImgSize batchSize]),'SSCB');
            if isempty(b) || (isnumeric(b) && all(b == 0))
                b = 0;
            end

            % Permute (H,W,In,Out) -> (H,W,Out,In) for dlconv
            Filter = permute(Filter, [1 2 4 3]);

            rImg = dlconv(inputImg,Filter,b, ...
                Stride=stride,DilationFactor=dilation,...
                    Padding=[padding(2) padding(1); padding(4) padding(3)]);
            r = reshape(extractdata(rImg),[],batchSize);
        else
            throw(CORAerror('CORA:notImplemented', ...
                'Non-dlconv execution not supported yet.'));
        end
    end

    function dW = convForWeigthsUpdate(obj,grad_out,input,options)
        % Compute dW via dlgradient: construct a pseudo-loss L = sum(Y.*dY)
        % where Y = dltranspconv(X, W, 0), then dW = dlgradient(L, W).

        [~,batchSize] = size(input);

        in_h = obj.inputSize(1);
        in_w = obj.inputSize(2);
        in_c = obj.inputSize(3);
        [out_h, out_w, out_c] = obj.aux_computeOutputSize();

        inputImg = dlarray(reshape(input,[in_h in_w in_c batchSize]),'SSCB');
        gradOutImg = dlarray(reshape(grad_out,[out_h out_w out_c batchSize]),'SSCB');
        W_dl = dlarray(obj.W);

        dW = extractdata(dlfeval(@gradient_closure, W_dl, inputImg, gradOutImg));

        function grad_W = gradient_closure(W_val, X_in, dY_in)
             % Forward pass with zero bias (bias is independent of dW)
             W_perm = permute(W_val, [1 2 4 3]);
             Y = dltranspconv(X_in, W_perm, 0, ...
                Stride=obj.stride, DilationFactor=obj.dilation, ...
                Cropping=[obj.padding(2) obj.padding(1); obj.padding(4) obj.padding(3)]);

             L = sum(Y .* dY_in, 'all');
             grad_W = dlgradient(L, W_val);
        end
    end


end

methods (Access = {?nnLayer, ?neuralNetwork})
    function bounds = evaluateInterval(obj, bounds, options)
        obj.checkInputSize()
        
        % Check input validity
        if any(bounds.inf > bounds.sup, 'all')
             throw(CORAerror('CORA:emptySet', ...
                 'Input interval to TransposedConv layer %s has inf > sup. Max violation: %e', ...
                 obj.name, max(bounds.inf - bounds.sup, [], 'all')));
        end

        % IBP (see Gowal et al. 2018)
        [mu,~] = obj.conv2d((bounds.sup + bounds.inf)/2,options, ...
            'sparseIdx');
        [r,~] = obj.conv2d((bounds.sup - bounds.inf)/2,options, ...
            'sparseIdx',abs(obj.W),[]);
            
        if any(r < 0, 'all')
             % CORAwarning('CORA:nnConvTranspose2DLayer:NegativeRadius', ...
             %    'Computed radius < 0 in layer %s. Min r: %e. Clamping to 0.', obj.name, min(r,[],'all'));
             r = max(r, 0);
        end

        l = mu - r;
        u = mu + r;
        bounds = interval(l,u);
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
        [nIn,~,~] = size(G);
        [nGrad,q,batchSize] = size(gG);

        if options.nn.interval_center
            cl = reshape(c(:,1,:),[nIn batchSize]);
            cu = reshape(c(:,2,:),[nIn batchSize]);
            gl = reshape(gc(:,1,:),[nGrad batchSize]);
            gu = reshape(gc(:,2,:),[nGrad batchSize]);

            % Backprop center interval
            [gl_, gu_] = backpropIntervalBatch(obj, cl, cu, gl, gu, options, updateWeights);
            gc_ = permute(cat(3,gl_,gu_),[1 3 2]);

            % Backprop generators
            [~,gG_] = obj.transconv2dZonotope(gl,gG,options,'sparseIdx',obj.W,[]);
        else
            if updateWeights
               [out_h,out_w,out_c] = obj.aux_computeOutputSize();
               biasUpdate = squeeze(sum(reshape(gc, ...
                    [out_h out_w out_c batchSize]),[1 2 4]));
               updateGrad(obj, 'b', biasUpdate, options);
            end
            % Backprop the center
            gc_ = backpropNumeric(obj, c, gc, options, updateWeights);

            % Backprop generators
            [~,gG_] = obj.transconv2dZonotope(gc,gG,options,'sparseIdx',obj.W,[]);
        end

        if updateWeights
             % Only using options.nn.train.zonotope_weight_update = 'sum'
             
             % G: (nIn, q, batchSize). (Small)
             G_flat = reshape(G, nIn, []); % (nIn, q*batchSize)
             % gG: (nOut, q, batchSize). (Large aka GradOut)
             gG_flat = reshape(gG, nGrad, []); % (nOut, q*batchSize)
             
             % Calculate weights update.
             % dW_gens = conv(gG_flat, G_flat_rotated).
             dW_gens = convForWeigthsUpdate(obj, gG_flat, G_flat, options);
             
             if ~options.nn.interval_center
                 % centerTerm = gc * c'.
                 dW_center = convForWeigthsUpdate(obj, gc, c, options);
                 weightsUpdate = dW_center + dW_gens;
             else
                 weightsUpdate = dW_gens;
             end
             
             updateGrad(obj, 'W', weightsUpdate, options);
        end
        
        % Set incoming center gradient.
        gc = gc_;
        gG = gG_;
    end
    % zonotope batch (for training / verification)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        obj.checkInputSize()

        if options.nn.interval_center
            [n,~,batchSize] = size(c);
            % Extract upper and lower bound.
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            
            % Sanitize bounds to prevent interval constructor crash
            if any(cl > cu, 'all')
                % CORAwarning('CORA:nnConvTranspose2DLayer:InvertedBounds', ...
                %    'Layer %s received inverted bounds (min cl-cu = %e). Swapping.', ...
                %    obj.name, min(cl-cu, [], 'all'));
                t = cl; 
                cl = min(t, cu); 
                cu = max(t, cu);
            end

            % Evaluate bounds.
            c = obj.evaluateInterval(interval(cl,cu),options);
            c = permute(cat(3,c.inf,c.sup),[1 3 2]);
            % Evaluate generators.
            c0 = zeros([n batchSize],'like',G);
            [~,G,Wff] = obj.conv2dZonotope(c0,G,options,'sparseIdx');
        else
            [c,G,Wff] = obj.conv2dZonotope(c,G,options,'sparseIdx');
        end

        % if options.nn.train.backprop
        %     obj.backprop.store.Wff = Wff;
        % end
    end
end

methods (Access = protected)
    function Wff = aux_conv2Mat(obj, varargin)
        % Compute matrix representation by transposing the equivalent
        % standard convolution matrix: Wff = M_conv'.

        [store, Filter, inImgSize, stride, padding, dilation] = ...
             setDefaultValues({'', obj.W, obj.inputSize, obj.stride, ...
                 obj.padding, obj.dilation}, varargin);

        [out_h, out_w, ~] = obj.aux_computeOutputSize(Filter, inImgSize, stride, padding, dilation);
        largeImgSize = [out_h, out_w, size(Filter,4)];

        % Build helper Conv2D with permuted weights (H,W,In,Out) -> (H,W,Out,In)
        W_helper = permute(Filter, [1 2 4 3]);
        b_helper = zeros(size(W_helper, 4), 1);
        helperLayer = nnConv2DLayer(W_helper, b_helper, padding, stride, dilation);
        helperLayer.inputSize = largeImgSize;

        % Transposed conv matrix is the transpose of the standard conv matrix
        Wff = helperLayer.aux_conv2Mat()';
    end

    function bias = aux_getPaddedBias(obj, varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

        % We need to expand bias to the Output Size (Large).
        [out_h, out_w, out_c] = obj.aux_computeOutputSize(Filter, ...
            inImgSize, stride, padding, dilation);

        if isempty(b)
            b = zeros(out_c,1,'like',Filter);
        elseif numel(b) == 1
            b = repmat(b(:),[out_c 1]);
        end

        % expand the bias vector to output size
        bias = repelem(b(:), out_h*out_w, 1);
    end
end

% Auxiliary functions -----------------------------------------------------

methods (Access = public)


    function outputSize = getOutputSize(obj, inImgSize)
        [out_h, out_w, out_c] = obj.aux_computeOutputSize(obj.W, inImgSize);
        outputSize = [out_h, out_w, out_c];
    end
end

end

% ------------------------------ END OF CODE ------------------------------
