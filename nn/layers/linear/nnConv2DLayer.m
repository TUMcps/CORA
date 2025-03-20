classdef nnConv2DLayer < nnLayer
% nnConv2DLayer - class for convolutional 2D layers
%
% Syntax:
%    obj = nnConv2DLayer(W, b, padding, stride, dilation, name)
%
% Inputs:
%    W - weight matrix (4-D single)
%        all the filters are stored in the weight-matrix:
%        (kernel_height, kernel_width, in_channels, num_filters)
%    b - bias column vector
%    padding - zero padding [left top right bottom]
%    stride - step size per dimension
%    dilation - step size per dimension
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] T. Gehr, et al. "AI2: Safety and Robustness Certification of
%        Neural Networks with Abstract Interpretation," 2018
%    [2] Practical Course SoSe '22 - Report Lukas Koller
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller, Tobias Ladner
% Written:       04-June-2022
% Last update:   01-December-2022 (combine with nnAvgPool2DLayer)
%                17-January-2023 (TL, optimizations)
%                13-December-2023 (LK, backpropagation + vectorized weight matrix generation)
% Last revision: 17-August-2022
%                02-December-2022 (clean up)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    W, b, stride, padding, dilation
end

methods
    % constructor
    function obj = nnConv2DLayer(varargin)
        
        % 1. parse input arguments: varargin -> vars
        narginchk(0,6)
        [W, b, padding, stride, dilation, name] = aux_parseInputArgs(varargin{:});
        
        % 2. check correctness of input arguments
        aux_checkInputArgs(W, b, padding, stride, dilation, name)

        % validate input
        [W, b, padding, stride, dilation, name] = setDefaultValues( ...
            {1, 0, [0 0 0 0], [1 1], [1 1], []}, varargin);
        inputArgsCheck({ ...
            {W, 'att', {'numeric', 'gpuArray'}}, ...
            {b, 'att', {'numeric', 'gpuArray'}}, ...
            {padding, 'att', 'numeric'}, ...
            {stride, 'att', 'numeric'}, ...
            {dilation, 'att', 'numeric'}, ...
        })

        % 3. call super class constructor
        obj@nnLayer(name)

        % 4. assign properties
        obj.W = W;
        obj.b = b;
        obj.padding = padding;
        obj.stride = stride;
        obj.dilation = dilation;

        % obj.backprop.store.Wff = [];
        % obj.backprop.store.WffIdx = [];
    end

    function [nin, nout] = getNumNeurons(obj)
        if isempty(obj.inputSize)
            nin = [];
            nout = [];
        else
            % we can only compute the number of neurons if the input
            % size was set.
            nin = obj.inputSize(1) * obj.inputSize(2) * obj.inputSize(3);
            outputSize = getOutputSize(obj, obj.inputSize);
            nout = prod(outputSize);
        end
    end

    % compute size of ouput feature map
    function outputSize = getOutputSize(obj, inImgSize)
        [out_h, out_w, out_c] = obj.aux_computeOutputSize(obj.W, inImgSize);
        outputSize = [out_h, out_w, out_c];
    end
end

% evaluate ------------------------------------------------------------

methods  (Access = {?nnLayer, ?neuralNetwork})
    % All dimensions (h,w,c) are flattened into a vector

    % numeric
    function r = evaluateNumeric(obj, input, options)
        obj.checkInputSize()
        [r,Wff] = obj.conv2d(input,'sparseIdx');

    end

    % numeric
    function bounds = evaluateInterval(obj, bounds, options)
        obj.checkInputSize()
        % IBP (see Gowal et al. 2018)
        [mu,~] = obj.conv2d((bounds.sup + bounds.inf)/2,'sparseIdx');
        % r = pagemtimes(abs(Wff),(bounds.sup - bounds.inf)/2);
        [r,~] = obj.conv2d((bounds.sup - bounds.inf)/2,'sparseIdx',abs(obj.W),[]);

        l = mu - r;
        u = mu + r;
        bounds = interval(l,u);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        obj.checkInputSize()

        % % compute weight and bias
        % Wff = obj.aux_conv2Mat();
        % bias = obj.aux_getPaddedBias();
        % 
        % % simulate using linear layer
        % linl = nnLinearLayer(Wff, bias);
        % Stmp = linl.evaluateSensitivity(S, x, options);

        [vK,vk,batchSize] = size(S);
        S = permute(S,[2 1 3]);
        S = reshape(S,[vk vK*batchSize]);
        S = obj.transconv2d(S,'sparseIdx',obj.W,[]);
        S = reshape(S,[size(x,1) vK batchSize]);
        S = permute(S,[2 1 3]);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        obj.checkInputSize()
        
        % compute weight and bias
        Wff = obj.aux_conv2Mat();
        bias = obj.aux_getPaddedBias();

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        [c, G, GI, E, id, id_, ind, ind_] = linl.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_, options);
    end

    % zonotope batch (for training)
    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        obj.checkInputSize()

        if options.nn.interval_center
            [n,~,batchSize] = size(G);
            % Extract upper and lower bound.
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            % Evaluate bounds.
            c = obj.evaluateInterval(interval(cl,cu));
            c = permute(cat(3,c.inf,c.sup),[1 3 2]);
            % Evaluate generators.
            c0 = zeros([n batchSize],'like',G);
            [~,G,Wff] = obj.conv2dZonotope(c0,G,'sparseIdx');
        else
            [c,G,Wff] = obj.conv2dZonotope(c,G,'sparseIdx');
        end

        % if options.nn.train.backprop
        %     obj.backprop.store.Wff = Wff;
        % end
    end

    % taylm
    function r = evaluateTaylm(obj, input, options)
        obj.checkInputSize()
    
        % compute weight and bias
        Wff = obj.aux_conv2Mat();
        bias = obj.aux_getPaddedBias();

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        r = linl.evaluateTaylm(obj, input, options);
    end

    % conZonotope
    function [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
        obj.checkInputSize()
        
        % compute weight and bias
        Wff = obj.aux_conv2Mat();
        bias = obj.aux_getPaddedBias();

        % simulate using linear layer
        linl = nnLinearLayer(Wff, bias);
        [c, G, C, d, l, u] = linl.evaluateConZonotope(obj, c, G, C, d, l, u, options);
    end

    % backprop ------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options)        
        % Compute weight update.
        dW = convForWeigthsUpdate(obj,grad_out,input);

        % Compute size of gradient.
        [out_h,out_w,out_c] = obj.aux_computeOutputSize();
        [~,batchSize] = size(input);

        % Compute bias update.
        db = squeeze(sum(reshape(grad_out, ...
            [out_h out_w out_c batchSize]),[1 2 4]));

        % Update weights and bias.
        updateGrad(obj, 'W', dW, options);
        updateGrad(obj, 'b', db, options);

        % The backproped gradient is computed by (full) convolving the 
        % outgoing gradient with the filters rotated by 180 degrees, which 
        % is the same as the transposed convolution.        
        grad_in = obj.transconv2d(grad_out,'sparseIdx',obj.W,[]);
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options)
        % See (Gowal et al. 2019)
        mu = (u + l)/2;
        r = (u - l)/2;

        % Compute weight update.
        dWmu = convForWeigthsUpdate(obj,gu + gl,mu);
        dWr = convForWeigthsUpdate(obj,gu - gl,r) .* sign(obj.W);

        % Compute size of gradient.
        [out_h,out_w,out_c] = obj.aux_computeOutputSize();
        [~,batchSize] = size(l);

        % Compute bias update.
        db = squeeze(sum(reshape(gl,[out_h out_w out_c batchSize]),[1 2 4]));

        % Update weights and bias.
        updateGrad(obj, 'W', dWmu + dWr, options);
        updateGrad(obj, 'b', db, options);

        % Use transposed convolutions to backprop gradients.
        dmu = obj.transconv2d((gu + gl)/2,'sparseIdx',obj.W,[]);
        dr = obj.transconv2d((gu - gl)/2,'sparseIdx',abs(obj.W),[]) ;
        gl = dmu - dr;
        gu = dmu + dr;
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options)
        in_h = obj.inputSize(1);
        in_w = obj.inputSize(2);
        in_c = obj.inputSize(3);
        
        % Compute size of gradient.
        [out_h,out_w,out_c] = obj.aux_computeOutputSize();
        [f_h, f_w] = obj.aux_getFilterSize();

        % padding [left,top,right,bottom]
        pad_l = obj.padding(1);
        pad_t = obj.padding(2);
        pad_r = obj.padding(3);
        pad_b = obj.padding(4);

        % Compute number of cropped rows and columns. A row or columns was
        % cropped during the forward propagation if the remaining number 
        % of input rows or columns does not fit the filter kernel.
        crop = ([in_h in_w] - [f_h f_w] + [pad_t + pad_b pad_l + pad_r]) ...
            - ([out_h out_w] - 1).*obj.stride;

        % gradient of the filters are computed by convolving the gradient 
        % with the input.
        [nIn,~,~] = size(G);
        [nGrad,q,batchSize] = size(gG);

        if options.nn.interval_center
            % Extract bounds.
            cl = reshape(c(:,1,:),[nIn batchSize]);
            cu = reshape(c(:,2,:),[nIn batchSize]);
            % Extract gradient for the bounds.
            gl = reshape(gc(:,1,:),[nGrad batchSize]);
            gu = reshape(gc(:,2,:),[nGrad batchSize]);
            [gl, gu] = backpropIntervalBatch(obj, cl, cu, gl, gu, options);
            gc = permute(cat(3,gl,gu),[1 3 2]);
            % % c = zeros(size(c.inf),'like',c.inf);
            % % gc = zeros(size(gc.inf),'like',gc.inf);
            % 
            % % Only using options.nn.train.zonotope_weight_update = 'sum'
            % % We move the generators to the batch. This is in order to do a
            % % convolution of the input with the outgoing gradient.
            % inputLin = reshape(G,nIn,[]);
            % inputLinPerm = reshape(permute(reshape(inputLin, ...
            %     [in_h in_w in_c q*batchSize]),[1 2 4 3]),[],in_c);
            % % Similarly, we move the generator of the outgoing gradient to the
            % % batch as well.
            % gradOutPermImg = permute(reshape(...
            %     reshape(gG,nGrad,[]),...
            %         [out_h out_w out_c q*batchSize]),[1 2 4 3]);
            % % To compute the weight update, the input is convolved with the
            % % outgoing gradient.
            % dW = obj.conv2d(inputLinPerm,'dWSparseIdx',gradOutPermImg,[],...
            %     [in_h in_w q*batchSize],obj.dilation,obj.padding,obj.stride);
            % [f_h, f_w] = obj.aux_getFilterSize();
            % weightsUpdate = permute(reshape(dW,[f_h f_w out_c in_c] + [crop 0 0]),[1 2 4 3]);
            % weightsUpdate = weightsUpdate(1:end-crop(1),1:end-crop(2),:,:);
            % % weightsUpdate = permute(reshape(dW,[f_h f_w out_c in_c]),[1 2 4 3]);
            % 
            % % % compute outer product of gradient and input zonotope
            % % genIds = obj.backprop.store.genIds;
            % % centerTerm = gc * c';
            % % gensTerm = 1/3*sum(pagemtimes(gG(:,genIds,:),'none', ...
            % %     G(:,genIds,:),'transpose'),3);
            % % 
            % % dW = reshape(accumarray( ...
            % %     reshape(obj.backprop.store.('sparseIdx'),[],1), ...
            % %     reshape(centerTerm + gensTerm,[],1)),[1 + f_h*f_w in_c out_c]);
            % % weightsUpdate = reshape(dW(2:end,:,:),[f_h f_w in_c out_c] + [crop 0 0]);
            % % weightsUpdate = weightsUpdate(1:end-crop(1),1:end-crop(2),:,:);
            % 
            % updateGrad(obj, 'W', weightsUpdate, options);
            % 
            % % The backproped gradient is computed by (full) convolving the 
            % % outgoing gradient with the filters rotated by 180 degrees, which 
            % % is the same as the transposed convolution.  
            c0 = zeros([nGrad batchSize],'like',gG);  
            [~,gG] = obj.transconv2dZonotope(c0,gG,'sparseIdx',obj.W,[]);
        else
            % Only using options.nn.train.zonotope_weight_update = 'sum'
            % We move the generators to the batch. This is in order to do a
            % convolution of the input with the outgoing gradient.
            inputLin = reshape(cat(2,permute(c,[1 3 2]),G),nIn,[]); % 1/3*
            inputLinPerm = reshape(permute(reshape(inputLin, ...
                [in_h in_w in_c (q+1)*batchSize]),[1 2 4 3]),[],in_c);
            % Similarly, we move the generator of the outgoing gradient to the
            % batch as well.
            gradOutPermImg = permute(reshape(...
                reshape(cat(2,permute(gc,[1 3 2]),gG),nGrad,[]),...
                    [out_h out_w out_c (q+1)*batchSize]),[1 2 4 3]);
            % To compute the weight update, the input is convolved with the
            % outgoing gradient.
            dW = obj.conv2d(inputLinPerm,'dWSparseIdx',gradOutPermImg,[],...
                [in_h in_w (q+1)*batchSize],obj.dilation,obj.padding,obj.stride);
            [f_h, f_w] = obj.aux_getFilterSize();
            weightsUpdate = permute(reshape(dW,[f_h f_w out_c in_c] + [crop 0 0]),[1 2 4 3]);
            weightsUpdate = weightsUpdate(1:end-crop(1),1:end-crop(2),:,:);
            % weightsUpdate = permute(reshape(dW,[f_h f_w out_c in_c]),[1 2 4 3]);
    
            % % compute outer product of gradient and input zonotope
            % genIds = obj.backprop.store.genIds;
            % centerTerm = gc * c';
            % gensTerm = 1/3*sum(pagemtimes(gG(:,genIds,:),'none', ...
            %     G(:,genIds,:),'transpose'),3);
            % 
            % dW = reshape(accumarray( ...
            %     reshape(obj.backprop.store.('sparseIdx'),[],1), ...
            %     reshape(centerTerm + gensTerm,[],1)),[1 + f_h*f_w in_c out_c]);
            % weightsUpdate = reshape(dW(2:end,:,:),[f_h f_w in_c out_c] + [crop 0 0]);
            % weightsUpdate = weightsUpdate(1:end-crop(1),1:end-crop(2),:,:);
    
            % Compute bias update.
            biasUpdate = squeeze(sum(reshape(gc, ...
                [out_h out_w out_c batchSize]),[1 2 4]));
            
            updateGrad(obj, 'W', weightsUpdate, options);
            updateGrad(obj, 'b', biasUpdate, options);
    
            % The backproped gradient is computed by (full) convolving the 
            % outgoing gradient with the filters rotated by 180 degrees, which 
            % is the same as the transposed convolution.    
            [gc,gG] = obj.transconv2dZonotope(gc,gG,'sparseIdx',obj.W,[]);
        end
    end
end

% Auxiliary functions ------------------------------------------------------

methods
    function layer = convert2nnLinearLayer(obj)
        % as convolutional layers are just fancy linear layers,
        % they can be converted to them easily.
        obj.checkInputSize()

        % compute weight and bias
        Wff = obj.aux_conv2Mat();
        bias = aux_getPaddedBias(obj);

        layer = nnLinearLayer(Wff, bias, sprintf("%s_linear", obj.name));
    end

    function names = getLearnableParamNames(obj)
        % list of learnable properties
        names = {'W', 'b'};
    end
end

methods (Access = protected)

    function [r,Wff] = conv2d(obj,input,varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

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
        rImg = dlconv(inputImg,Filter,b, ...
            Stride=stride,DilationFactor=dilation,...
                Padding=[pad_t pad_l; pad_b pad_r]);
        r = reshape(extractdata(rImg),[],batchSize);

        Wff = [];

        % % compute weight and bias
        % bias = obj.aux_getPaddedBias(varargin{:});
        % 
        % if length(varargin) >= 3
        %     varargin = varargin([1:2,4:end]);
        % end
        % 
        % % compute weight matrix
        % Wff = obj.aux_conv2Mat(varargin{:});
        % r = Wff * input + bias;
        % % r = Wff * sparse(double(input));
        % % r = single(full(r)) + bias;

    end

    function [c,G,Wff] = conv2dZonotope(obj,c,G,varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

        % Put generators into batch and do regular convolution.
        [n,q,batchSize] = size(G);
        inputLin = reshape(cat(2,permute(c,[1 3 2]),G),n,(q+1)*batchSize);

        [rLin,Wff] = obj.conv2d(inputLin,varargin{:});
        r = reshape(rLin,[],q+1,batchSize);

        c = reshape(r(:,1,:),[size(r,1) batchSize]);
        G = r(:,2:end,:);

        % % compute weight and bias
        % bias = obj.aux_getPaddedBias(varargin{:});
        % 
        % if length(varargin) >= 3
        %     varargin = varargin([1:2,4:end]);
        % end
        % 
        % % compute weight matrix
        % Wff = obj.aux_conv2Mat(varargin{:});
        % c = Wff*c + bias;
        % G = pagemtimes(Wff,G);
    end

    function r = transconv2d(obj,input,varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ... 
                obj.stride, obj.padding, obj.dilation}, varargin);

        [~,batchSize] = size(input);

        % padding [left,top,right,bottom]
        pad_l = padding(1);
        pad_t = padding(2);
        pad_r = padding(3);
        pad_b = padding(4);

        % Compute size of gradient.
        [out_h,out_w,out_c] = obj.aux_computeOutputSize(Filter,inImgSize,stride,...
            padding, dilation);

        inputImg = dlarray(reshape(input,[out_h out_w out_c batchSize]),'SSCB');
        if isempty(b)
            b = 0;
        end
        rImg = dltranspconv(inputImg,Filter,b, ...
            Stride=stride,DilationFactor=dilation,...
                Cropping=[pad_t pad_l; pad_b pad_r]);

        in_w = inImgSize(2);
        in_c = inImgSize(3);
        % compute cropped pixels and insert zeros bottom and right
        crop = obj.computeCrop();
        rImg = [rImg zeros([size(rImg,1) crop(2) in_c batchSize],'like',rImg);
            zeros([crop(1) in_w in_c batchSize],'like',rImg)];
        r = reshape(extractdata(rImg),[],batchSize);

        % if length(varargin) >= 3
        %     varargin = varargin([1,4:end]);
        % end
        % 
        % % compute weight matrix
        % Wfft = obj.aux_conv2Mat(varargin{:})';
        % r = Wfft * input;
        % % r = Wfft * sparse(double(input));
        % % r = single(full(r)); % + bias;
    end

    function [c,G] = transconv2dZonotope(obj,c,G,varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.inputSize, obj.W, obj.b, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

        % Put generators into batch and do regular convolution.
        [n,q,batchSize] = size(G);
        inputLin = reshape(cat(2,permute(c,[1 3 2]),G),n,[]);

        rLin = obj.transconv2d(inputLin,varargin{:});
        r = reshape(rLin,[],q+1,batchSize);

        c = squeeze(r(:,1,:));
        G = r(:,2:end,:);

        % if length(varargin) >= 3
        %     varargin = varargin([1,4:end]);
        % end
        % 
        % Wfft = obj.aux_conv2Mat(varargin{:})';
        % c = Wfft*c;
        % G = pagemtimes(Wfft,G);
    end

    function crop = computeCrop(obj)
        in_h = obj.inputSize(1);
        in_w = obj.inputSize(2);
        
        % Compute size of gradient.
        [out_h,out_w,~] = obj.aux_computeOutputSize();

        % padding [left,top,right,bottom]
        pad_l = obj.padding(1);
        pad_t = obj.padding(2);
        pad_r = obj.padding(3);
        pad_b = obj.padding(4);

        % Compute number of cropped rows and columns. A row or columns was
        % cropped during the forward propagation if the remaining number 
        % of input rows or columns does not fit the filter kernel.
        [f_h, f_w] = obj.aux_getFilterSize();
        crop = ([in_h in_w] - [f_h f_w] + [pad_t + pad_b pad_l + pad_r]) ...
            - ([out_h out_w] - 1).*obj.stride;
    end

    function dW = convForWeigthsUpdate(obj,grad_out,input)
        [~,batchSize] = size(input);
        in_h = obj.inputSize(1);
        in_w = obj.inputSize(2);
        in_c = obj.inputSize(3);
        
        % Compute size of gradient.
        [out_h,out_w,out_c] = obj.aux_computeOutputSize();

        % padding [left,top,right,bottom]
        % pad_l = obj.padding(1);
        % pad_t = obj.padding(2);
        % pad_r = obj.padding(3);
        % pad_b = obj.padding(4);

        % Compute number of cropped rows and columns. A row or columns was
        % cropped during the forward propagation if the remaining number 
        % of input rows or columns does not fit the filter kernel.
        % [f_h, f_w] = obj.aux_getFilterSize();
        % crop = ([in_h in_w] - [f_h f_w] + [pad_t + pad_b pad_l + pad_r]) ...
        %     - ([out_h out_w] - 1).*obj.stride;
        crop = obj.computeCrop();

        % Alternative implementation: compute weights update with
        % 'accumarray'.
        % WffIdx = obj.backprop.store.WffIdx;
        % dW = reshape(accumarray( ...
        %     reshape(WffIdx,[],1), ...
        %     reshape(grad_out * input',[],1)),[1 + f_h*f_w in_c out_c]);
        % weightsUpdate0 = reshape(dW(2:end,:,:),[f_h f_w in_c out_c]);

        % To compute the weight update, the input is convolved with the
        % outgoing gradient.
        inputPerm = reshape(permute(reshape(input, ...
            [in_h in_w in_c batchSize]),[1 2 4 3]),[],in_c);
        gradOutPermImg = permute(reshape(grad_out, ...
            [out_h out_w out_c batchSize]),[1 2 4 3]);
        % Compute weight update.
        dW = obj.conv2d(inputPerm,'dWSparseIdx',gradOutPermImg,[],...
            [in_h in_w batchSize],obj.dilation,obj.padding,obj.stride);
        % Reshape and crop to size.
        [f_h, f_w] = obj.aux_getFilterSize();
        dW = permute(reshape(dW,[f_h f_w out_c in_c] + [crop 0 0]),[1 2 4 3]);
        dW = dW(1:end-crop(1),1:end-crop(2),:,:);
    end

    function [f_h, f_w] = aux_getFilterSize(obj, varargin)
        [Filter, dilation] = setDefaultValues({obj.W, obj.dilation}, ...
            varargin);
        % return the size of the filter kernels.
        k_h = size(Filter, 1);
        k_w = size(Filter, 2);
        d_h = dilation(1);
        d_w = dilation(2);

        f_h = k_h + (k_h - 1) * (d_h - 1);
        f_w = k_w + (k_w - 1) * (d_w - 1);
    end

    function [out_h, out_w, out_c] = aux_computeOutputSize(obj, varargin)
        [Filter, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({obj.W, obj.inputSize, obj.stride, ...
                obj.padding, obj.dilation}, varargin);

        in_h = inImgSize(1);
        in_w = inImgSize(2);
        [f_h, f_w] = aux_getFilterSize(obj, Filter, dilation);
        % padding [left,top,right,bottom]
        pad_l = padding(1);
        pad_t = padding(2);
        pad_r = padding(3);
        pad_b = padding(4);
        % stride
        stride_h = stride(1);
        stride_w = stride(2);

        out_h = floor((in_h - f_h + pad_t + pad_b)/stride_h) + 1;
        out_w = floor((in_w - f_w + pad_l + pad_r)/stride_w) + 1;
        out_c = size(Filter, 4);
    end

    % Compute an index-matrix to express convolutions as matrix-vector multiplication.
    function [WffIdx,WffSparseIdx] = aux_computeWeightMatIdx(obj, varargin)
        [store, Filter, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.inputSize, obj.stride, ...
                obj.padding, obj.dilation}, varargin);
        
        if ~isempty(store) && isfield(obj.backprop.store,store)
            % Return stored index matrix.
            WffIdx = obj.backprop.store.(store);

            WffSparseIdx = {}; % obj.backprop.store.(store);
            return
        end

        % Here we compute the index-matrix to express a convolution as a 
        % matrix-vector multiplication.

        % The construction of the index-matrix is vectorized. The
        % implementation is based on the MATLAB toeplitz function. The
        % entire weight matrix can be decomposed into weight matrices one
        % for each pair of input-output channel. Each smaller weight matrix
        % has its filter weights at the exact same positions. Thus, each
        % filter matrix is constructed using the same index matrix, which
        % indexes the filter kernel.
        % An index matrix for the entire weight matrix is obtained by
        % repeating the small index matrix and adding a corresponding
        % offset for each input-output-channel pair.
        % The small index matrix is constructed by first flattening the
        % indices (1:k_h*k_w) into a row-vector. Zeros are pre- and 
        % appended s.t. the entries (end-in_h*in_w:end) of the resulting 
        % row-vector form the first line of the index matrix and the 
        % entries (1:in_h*in_w) correspond to the last row of the index 
        % matrix. The idea now is to use MATLAB's vector-addition of a 
        % column-vector and a row-vector to create the index matrix. The
        % column-vector determines the shift of the filter to the next row,
        % while the row-vector is just the flattened & zero-padded indices.

        % Finally, we obtain a matrix containing indices for the
        % filter-kernels and zeros. We shift the indices by 1 and prepend a
        % single 0 to the filter-kernels.

        % Obtain input image size.
        in_h = inImgSize(1);
        in_w = inImgSize(2);

        % Compute output size.
        [out_h, out_w, ~] = obj.aux_computeOutputSize(Filter, ...
            inImgSize, stride, padding, dilation);

        % Get padding [left,top,right,bottom].
        pad_l = padding(1);
        pad_t = padding(2);
        pad_r = padding(3);
        pad_b = padding(4);
        % Pad image with zeros.
        in_h_pad = in_h + pad_t + pad_b;
        in_w_pad = in_w + pad_l + pad_r;

        % Get stride values.
        s_h = stride(1);
        s_w = stride(2);

        n = out_h*out_w;
        m = in_h_pad*in_w_pad;

        % Get filter kernel size.
        [k_h, k_w, in_c, out_c] = size(Filter);

        % Compute an index-matrix for an individual filter.
        filterIdx = cast(1:k_h*k_w,'like',Filter);
        filterIdx = reshape(filterIdx,k_h,k_w);
        if any(dilation ~= 1)
            % Dilate filter indices, by padding them with zeros and permuting
            % the rows and columns. Get dilation factors.
            d_h = dilation(1);
            d_w = dilation(2);
            % Compute dilated filter size.
            f_h = k_h + (k_h - 1) * (d_h - 1);
            f_w = k_w + (k_w - 1) * (d_w - 1);
            % Append zero-rows and -columns.
            padFilterIdx = eye(f_h,k_h)*filterIdx*eye(k_w,f_w);
            % Permute the rows and columns s.t. appended zero-rows and -columns 
            % are inbetween the filter entries.
            idx_h = reshape([1:k_h-1; ...
                reshape(k_h+1:(k_h-1)*d_h+1,k_h-1,d_h-1)'],[],1);
            idx_w = reshape([1:k_w-1; ...
                reshape(k_w+1:(k_w-1)*d_w+1,k_w-1,d_w-1)'],[],1);
            filterIdx = padFilterIdx([idx_h; k_h],[idx_w; k_w]);
        else
            % Filter is not dilated; filter size is equal to the kernel
            % size.
            f_h = k_h;
            f_w = k_w;
        end

        % We turn the filter kernel into a vector, padded with 0's.
        linFilterIdx = [filterIdx ...
            zeros(f_h,in_w_pad-f_w,'like',Filter); 
            zeros(in_h_pad-f_h,in_w_pad,'like',Filter)];
        % Compute number of zeros needed for padding before the filter.
        nzspad = in_h_pad*(in_w_pad-f_w+1)-f_h;
        linFilterIdx = [
            zeros(nzspad,1,'like',Filter); 
            linFilterIdx(:)];
        % We compute indices that index the vectorized filter to obtain the
        % filter-matrix. We compute the shifts for each rows; thereby, we
        % have to adjust for row breaks and stride.
        % Adjust for row breaks.
        rowShift = cast(out_w-1:-1:0,'like',Filter);
        rowShift = reshape(repmat(rowShift,out_h,1),[],1); 
        % Adjust for horizontal and vertical stride.
        rowShift = ((s_w-1)*in_h_pad + f_w-1)*rowShift + ... % horizontal stride
            (s_h-1)*(repmat(out_h-1:-1:0,1,out_w)' + (out_w-1)*rowShift); % vertical stride
        % Adjust for cutoff rows and columns. Happens when filter does not
        % fit for the last rows or columns.
        rowShift = rowShift + ...
            in_h_pad*mod(in_h_pad - f_h,s_h) + ... % adjust for cutoff rows
            mod(in_w_pad - f_w,s_w)*reshape(repmat(out_w:-1:1,out_h,1),[],1); % adjust for cutoff columns
        % Compute the index matrix.
        ij = (cast((n:-1:1)','like',Filter) + rowShift) + cast(0:m-1,'like',Filter);
        % ij = (cast((nzspad:-1:(nzspad-n+1))','like',Filter) + rowShift) + cast(0:m-1,'like',Filter);
        % Compute the filter-index matrix. Has the shape of the filter
        % matrix and indexes into a filter kernel.
        WfIdx = reshape(linFilterIdx(ij),[n m]);
        % We remove padded area.
        isPad = false(in_h_pad, in_w_pad);
        isPad(1:pad_t, :) = true;
        isPad(end+1-pad_b:end, :) = true;
        isPad(:, 1:pad_l) = true;
        isPad(:, end+1-pad_r:end) = true;
        isPad = reshape(isPad, 1, []);
        WfIdx(:, isPad) = [];
        % We add 1 to allow for indeing a zero value which is prepended to 
        % each filter.
        WfIdx = WfIdx + 1;

        % Assemble individual filter-index matrices to larger index matrix.
        WffIdx = repmat(WfIdx,out_c,in_c);
        % We have to shift the indices for each filter by the number of 
        % entries in each filter.
        indShift = (k_h*k_w + 1)*cast(0:in_c*out_c-1,'like',Filter);
        indShift = repelem(reshape(indShift,in_c,out_c)',out_h*out_w,in_h*in_w);
        WffIdx = WffIdx + indShift;

        % % compute sparse indices
        % WfIdxSparse = sparse(double(WfIdx - 1));
        % [i,j,s] = find(WfIdxSparse);
        % % number of nonzero indices
        % numNonzero = nnz(WfIdxSparse);
        % 
        % % We have to shift the indices for each filter by the number of 
        % % entries in each filter.
        % 
        % colShift = repmat(cast(0:out_c-1,'like',Filter),1,in_c)';
        % colShift = repelem(colShift,numNonzero,1);
        % i = repmat(i,in_c*out_c,1) + out_h*out_w*colShift;
        % 
        % rowShift = repelem(cast(0:in_c-1,'like',Filter),1,out_c)';
        % rowShift = repelem(rowShift,numNonzero,1);
        % j = repmat(j,in_c*out_c,1) + in_h*in_w*rowShift;
        % 
        % indShift = cast(0:in_c*out_c-1,'like',Filter)';
        % indShift = repelem(indShift,numNonzero,1);
        % s = repmat(s,in_c*out_c,1) + k_h*k_w*indShift;
        % 
        % n = size(WfIdx,1)*out_c;
        % m = size(WfIdx,2)*in_c;
        % 
        % WffSparseIdx = {i,j,s,n,m};
        WffSparseIdx = {};

        if ~isempty(store)
            % Store index matrix
            % obj.backprop.store.WffIdx = WffIdx;

            obj.backprop.store.(store) = WffIdx; % WffSparseIdx;
        end
    end

    % Compute weight matrix to express convolutions as matrix-vector multiplication.
    function Wff = aux_conv2Mat(obj, varargin)
        [store, Filter, inImgSize, stride, padding, dilation] = ... 
            setDefaultValues({'', obj.W, obj.inputSize, obj.stride, ... 
                obj.padding, obj.dilation}, varargin);

        % Get filter kernel size.
        [~, ~, in_c, out_c] = size(Filter);

        % WffIdx = obj.aux_computeWeightMatIdx(varargin{:});
        [WffIdx,WffSparseIdx] = obj.aux_computeWeightMatIdx(varargin{:});
        % Vectorize all filter kernels and prepend 0.
        linFfilter = [zeros(1,in_c,out_c,'like',Filter); 
            reshape(Filter,[],in_c,out_c)];
        linFfilter = linFfilter(:);
        % Construct weight matrix.
        Wff = reshape(linFfilter(WffIdx),size(WffIdx));

        % Alternative implementation with sparse matrices (also see 
        % aux_computeWeightMatIdx).
        % WffSparseIdx{3} = double(Filter(WffSparseIdx{3}));
        % WffSparse = sparse(WffSparseIdx{:});
        % Wff = single(full(WffSparse));
    end

    % Pad bias s.t. convolution can be simulated by linear layer.
    function bias = aux_getPaddedBias(obj, varargin)
        [store, Filter, b, inImgSize, stride, padding, dilation] = ...
            setDefaultValues({'', obj.W, obj.b, obj.inputSize, ...
                obj.stride, obj.padding, obj.dilation}, varargin);

        % We turn the bias into a vector.

        % compute output size
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

end


% Auxiliary functions -----------------------------------------------------

function [W, b, padding, stride, dilation, name] = aux_parseInputArgs(varargin)
    % validate input
    [W, b, padding, stride, dilation, name] = setDefaultValues( ...
        {1, 0, [0 0 0 0], [1 1], [1 1], []}, varargin);
end

function aux_checkInputArgs(W, b, padding, stride, dilation, name)
    % check input types
    inputArgsCheck({ ...
        {W, 'att', 'numeric'}, ...
        {b, 'att', 'numeric'}, ...
        {padding, 'att', 'numeric'}, ...
        {stride, 'att', 'numeric'}, ...
        {dilation, 'att', 'numeric'}, ...
        % name is checked in nnLayer
    })

    % check dimensions
    if length(size(W)) > 4
        throw(CORAerror('CORA:wrongInputInConstructor','Weight matrix has wrong dimensions.'))
    end
    if size(W,4) ~= length(b)
        throw(CORAerror("CORA:wrongInputInConstructor",'Weight matrix and bias dimensions don''t match.'))
    end
    if length(padding) ~= 4
        throw(CORAerror('CORA:wrongInputInConstructor', 'Padding must be an array with 4 entries: [left top right bottom]'))
    end
    if length(stride) ~= 2
        throw(CORAerror('CORA:wrongInputInConstructor', 'Stride must be an array with 2 entries: [dim1, dim2]'))
    end
    if length(dilation) ~= 2
        throw(CORAerror('CORA:wrongInputInConstructor', 'Dilation must be an array with 2 entries: [dim1, dim2]'))
    end

    
end

% ------------------------------ END OF CODE ------------------------------
