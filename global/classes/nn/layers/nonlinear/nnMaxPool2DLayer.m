classdef nnMaxPool2DLayer < nnLayer
% nnMaxPool2DLayer - class for max pooling 2D layers
%
% Syntax:
%    obj = nnMaxPool2DLayer(poolSize, stride, name)
%
% Inputs:
%    poolSize - size of pooling area, column vector
%    stride - step size, column vector
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
% See also: NeuralNetwork

% Authors:       Lukas Koller, Tobias Ladner
% Written:       05-June-2022
% Last update:   02-December-2022 (better permutation computation)
% Last revision: 17-August-2022
%                02-December-2022 (clean up)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
end

properties
    poolSize, stride %, padding
end

methods
    % constructor
    function obj = nnMaxPool2DLayer(poolSize, stride, name)
        % TODO padding, dilation

        if nargin < 2
            stride = poolSize;
        end
        if nargin < 3
            name = [];
        end
        % call super class constructor
        obj@nnLayer(name)

        obj.poolSize = poolSize;
        obj.stride = stride;
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

    % compute output size
    function outputSize = getOutputSize(obj, imgSize)
        in_h = imgSize(1);
        in_w = imgSize(2);
        pool_h = obj.poolSize(1);
        pool_w = obj.poolSize(2);
        % if pool width or height do not divide image width or height
        % the remaining pixels are ignored
        out_h = floor((in_h - pool_h)/obj.stride(1)) + 1;
        out_w = floor((in_w - pool_w)/obj.stride(2)) + 1;
        out_c = imgSize(3); % preserve number of channels
        outputSize = [out_h, out_w, out_c];
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    % input: column vector, all dimensions (h,w,c) are flattened into a vector

    % numeric
    function r = evaluateNumeric(obj, input, evParams)
        num_samples = size(input, 2);

        % move adjacent pixels next to each other
        id_mpp = aux_computePermutationMatrix(obj);

        % rearrange
        input = input(id_mpp, :);
        input = reshape(input, prod(obj.poolSize), [], num_samples);

        % compute max
        r = max(input, [], 1);

        % reshape to match input
        r = reshape(r, [], num_samples);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, evParams)
        % compute id with maximum value
        function [dmax, maxId] = getMaxId(ids)
            [~, argmax] = max(x(ids));
            maxId = ids(argmax);
            dmax = 0;
        end

        % compute projection matrix
        [~, P] = aux_computeLinearProjectionMatrix(obj, @getMaxId);
        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        bias = zeros(prod(outputSize), 1);

        % simulate using linear layer
        linl = nnLinearLayer(P, bias);
        S = linl.evaluateSensitivity(S, x);
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        checkInputSize(obj)

        switch evParams.maxpool_type
            case 'project'
                % Projection maximum pooling approach: see [2]
                % 1. Determine dimensions where pooling is applied
                % 2. Choose 'most likely' max dimension
                % 3. Project onto this dimension
                % 4. Add approximation error to polynomial zonotope
                [c, G, GI, E, id, id_, ind, ind_] = aux_evaluatePolyZonotopeProject(obj, c, G, GI, E, id, id_, ind, ind_, evParams);
            case 'regression'
                % Regression maximum pooling approach: see [2]
                % 1. find a polynomial to approximate maximum pooling.
                %    For each pooled area we find a polynomial.
                % 2. evaluate polynomial with polynomial zonotope
                % 3. bound approximation error
                % 4. add approximation error to polynomial zonotope
                [c, G, GI, E, id, id_, ind, ind_] = aux_evaluatePolyZonotopeRegression(obj, c, G, GI, E, id, id_, ind, ind_, evParams);
            otherwise
                throw(CORAerror('CORA:wrongFieldValue', 'evParams.maxpool_type', {'project', 'regression'}))
        end
    end
end

% internal functions ------------------------------------------------------

methods (Access = private)
    % Compute permutation-matrix to permute an input vector s.t.
    % pool-ed elements are adjacent. Only computes permutation matrix
    % for a single channel. This matrix also removes row and columns
    % s.t. the pooling size divides the image size. see [2]
    function Wmp = aux_computePermutationMatrixForChannel(obj)
        % see [1]
        % Only compute the permutation-matrix for a single channel.
        img_h = obj.inputSize(1);
        img_w = obj.inputSize(2);

        pool_h = obj.poolSize(1);
        pool_w = obj.poolSize(2);

        num_pools_h = floor(img_h/pool_h);
        num_pools_w = floor(img_w/pool_w);

        pools_h = eye(num_pools_h);
        unitvector = @(i, n) (1:n == i) * 1;
        pools_hw = cell(1, pool_w);
        for j = 1:pool_w
            pools_hj = kron(pools_h, unitvector(j, pool_w)');
            entries_hj = kron(pools_hj, eye(pool_h));
            entries_hj = [; ...
                entries_hj, ...
                zeros(size(entries_hj, 1), img_h-(pool_h * num_pools_h)); ...
                ];

            pools_hw{j} = entries_hj;
        end

        pools_hw = [pools_hw{:}];
        Wmp = kron(eye(num_pools_w), pools_hw);
        Wmp = [Wmp, zeros(size(Wmp, 1), img_h*(img_w - (pool_w * num_pools_w)))];

    end

    % Compute permutation-matrix for entinre input vector. see [2]
    function id_mpp = aux_computePermutationMatrix(obj)
        % compute permutation matrix for single channel
        Wmp = aux_computePermutationMatrixForChannel(obj);

        % construct permutation matrix for entire input
        c_in = obj.inputSize(3);
        id_mpp = 1:prod(obj.inputSize);
        id_mpp = reshape(id_mpp, [], c_in);
        id_mpp = Wmp * id_mpp;
        id_mpp = reshape(id_mpp, [], 1);
    end

    % Compute projection matrix to simulate max-pool with a linear
    % layer and error per dimension. Returns the projection matrix and
    % a vector with the approximation error for each output dimension.
    % see [2]
    function [d, P] = aux_computeLinearProjectionMatrix(obj, getMaxId)
        % [dmax,maxId] = getMax(ids): function that computes id with
        % maximum value and potential error
        img_h = obj.inputSize(1);
        img_w = obj.inputSize(2);
        c_in = obj.inputSize(3);

        pool_h = obj.poolSize(1);
        pool_w = obj.poolSize(2);
        num_pools_h = floor(img_h/pool_h); % number of pools in fst dim
        num_pools_w = floor(img_w/pool_w); % number of pools in snd dim

        % compute permutated indices to rearrange pixels in the same pool
        % next to each other
        id_mpp = aux_computePermutationMatrix(obj);

        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        img_h_out = outputSize(1);
        img_w_out = outputSize(2);
        c_out = outputSize(3);
        % init projection matrix
        P = zeros(img_h_out*img_w_out*c_out, img_h*img_w*c_in);
        % init error vector
        d = zeros(img_h_out*img_w_out*c_out, 1);
        for k = 1:c_out % number of out-channels
            for i = 1:(num_pools_h * num_pools_w) % iterate over pool areas
                % compute the start index of the pool-ed area
                poolStart = (k - 1) * num_pools_h * pool_h * num_pools_w * pool_w + (i - 1) * pool_h * pool_w + 1;
                % find max element
                poolIds = poolStart:(poolStart + pool_h * pool_w - 1);
                [dmax, maxId] = getMaxId(id_mpp(poolIds));
                % set max element in output
                outId = (k - 1) * img_h_out * img_w_out + i;
                P(outId, maxId) = 1;
                d(outId) = dmax;
            end
        end
    end

    % evParams.maxpool_type = 'project'
    function [c, G, GI, E, id, id_, ind, ind_] = aux_evaluatePolyZonotopeProject(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        % input: column vector, all dimensions (h,w,c) are flattened into a vector
        checkInputSize(obj)

        % compute the interval for each dimension
        pZ = polyZonotope(c, G, GI, E, id);
        I = interval(pZ); % ,'split');

        % compute id with maximum value and error
        function [dmax, maxId] = aux_getMaxId(ids)
            % compute max infimum
            maxInf = max(infimum(I(ids)));
            % only keep ids with supremum greater than maxInf
            ids = ids(supremum(I(ids)) >= maxInf);
            ivals = I(ids);

            A = (ones(length(ids)) - eye(length(ids))) .* (repmat(supremum(ivals), 1, length(ids)).') - infimum(ivals);
            A(1:length(ids)+1:end) = -Inf;
            [dmax, maxId] = min(max(A.'));
            dmax = max(0, dmax);
            maxId = ids(maxId);
        end

        % compute projection matrix
        [d, P] = aux_computeLinearProjectionMatrix(obj, @aux_getMaxId);
        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        bias = zeros(prod(outputSize), 1);

        % simulate using linear layer
        linl = nnLinearLayer(P, bias);
        [c, G, GI, E, id, id_, ind, ind_] = linl.evaluatePolyZonotope(c, G, GI, E, id, id_, ind, ind_);

        c = c + 0.5 * d;

        % add error d
        Gd = diag(0.5*d);
        % Gd = diag(d);
        if evParams.add_approx_error_to_GI
            GI = [GI, Gd];
        else
            G = [G, Gd];
            E = blkdiag(E, eye(size(Gd, 2)));
            id = [id; 1 + (1:size(Gd, 2))' * id_];
            id_ = max(id);
        end
    end

    % evParams.maxpool_type = 'regression'
    function [c, G, GI, E, id, id_, ind, ind_] = aux_evaluatePolyZonotopeRegression(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
        % set polynom order
        order = 1; % only works for order=1!

        checkInputSize(obj)
        img_h = obj.inputSize(1);
        img_w = obj.inputSize(2);

        pool_h = obj.poolSize(1);
        pool_w = obj.poolSize(2);
        num_pools_h = floor(img_h/pool_h); % number of pools in fst dim
        num_pools_w = floor(img_w/pool_w); % number of pools in snd dim


        % compute permutated indices to rearrange pixels in the same pool
        % next to each other
        id_mpp = aux_computePermutationMatrix(obj);

        % compute interval for each dimension
        pZ = polyZonotope(c, G, GI, E, id);
        I = interval(pZ);

        % compute output size
        outputSize = getOutputSize(obj, obj.inputSize);
        img_h_out = outputSize(1);
        img_w_out = outputSize(2);
        c_out = outputSize(3);

        % init output
        c_ = zeros(img_h_out*img_w_out*c_out, 1);
        G_ = zeros(img_h_out*img_w_out*c_out, size(G, 2));
        GI_ = zeros(img_h_out*img_w_out*c_out, size(GI, 2));
        d = zeros(img_h_out*img_w_out*c_out, 1);

        for k = 1:c_out % number of out-channels
            for i = 1:(num_pools_h * num_pools_w) % iterate over pool areas
                % compute the start index of the pool-ed area
                poolStart = (k - 1) * num_pools_h * pool_h * num_pools_w * pool_w + (i - 1) * pool_h * pool_w + 1;
                % compute ids of pool-ed area in input vector
                poolIds = id_mpp(poolStart:(poolStart + pool_h * pool_w - 1));

                % 0. exclude dimensions that can never become maximum.

                % compute max infimum
                maxInf = max(infimum(I(poolIds)));
                % only keep ids with supremum greater than maxInf
                poolIds = poolIds(supremum(I(poolIds)) >= maxInf);
                ivals = I(poolIds);
                % id in output
                outId = (k - 1) * img_h_out * img_w_out + i;

                % when there is only one dimension remaining, take this dimension
                if length(poolIds) == 1
                    c_(outId) = c(poolIds);
                    G_(outId, :) = G(poolIds, :);
                    if ~isempty(GI)
                        GI_(outId, :) = GI(poolIds, :);
                    end
                    continue;
                end

                % 1. Find polynomial to approximate maximum of the
                % pooled area.

                % sample points
                samples = 100;
                x = zeros(samples, length(poolIds));
                for j = 1:length(poolIds)
                    x(:, j) = linspace(infimum(ivals(j)), supremum(ivals(j)), samples);
                end
                % compute max accross the samples
                y = max(x, [], 2);
                % set up X matrix. X*coeffs = y
                X = zeros(samples, (order + 1)*length(poolIds));
                for j = 1:length(poolIds)
                    X(:, (j - 1)*(order + 1)+1:j*(order + 1)) = x(:, j).^(0:order);
                end
                % find coefficients
                coeffs = pinv(X) * y;
                coeffs = reshape(coeffs, [order + 1, length(poolIds)]);
                coeffs = coeffs';

                % 2. evaluate the polynomial with polynomial zonotope

                % evaluate polynomial on poly zonotope
                for j = 1:length(poolIds)
                    cj = c(poolIds(j));
                    Gj = G(poolIds(j), :);
                    if ~isempty(GI)
                        GIj = GI(poolIds(j), :);
                    end

                    c_(outId) = c_(outId) + coeffs(j, 1) + coeffs(j, 2) * cj;
                    G_(outId, :) = G_(outId, :) + coeffs(j, 2) * Gj;
                    if ~isempty(GI)
                        GI_(outId, :) = GI_(outId, :) + coeffs(j, 2) * GIj;
                    end
                end

                % 3. bound approximation error

                % find max difference
                num_samples = 10000; % in total sample ~10000 points
                % divide number of samples accross dimensions
                dnom = sum(supremum(ivals)-infimum(ivals));
                n = zeros(length(poolIds), 1);
                for j = 1:length(poolIds)
                    sup = supremum(ivals(j));
                    inf = infimum(ivals(j));
                    if dnom == 0
                        n(j) = 1;
                    else
                        e = log10(num_samples) * ((sup - inf) / dnom);
                        n(j) = ceil(10^e);
                    end
                end
                % compute tolerances for each dimension (max derBounds are [0,1])
                der1 = interval(0, 1);
                der = zeros(length(poolIds), 1);
                for j = 1:length(poolIds)
                    der2 = interval(coeffs(j, 2), coeffs(j, 2)); % only linear function; no need to call nnHelper.getDerInterval(coeffs(j,:),inf,sup);
                    %der2 = -der2;
                    der(j) = supremum(abs(der1-der2));
                end
                dx = (supremum(ivals) - infimum(ivals)) ./ n; % stepsize per dimension
                tol = der .* dx; % maximum change in between samples per dimension
                % sample points uniformly in all dimensions
                x = zeros(length(poolIds), prod(n));
                for j = 1:length(poolIds)
                    sup = supremum(ivals(j));
                    inf = infimum(ivals(j));
                    xj = linspace(inf, sup, n(j));
                    x(j, :) = reshape(repmat(xj, prod(n(1:j-1)), prod(n(j+1:end))), [], 1);
                end

                % flip coeffs for polyVal
                coeffs = fliplr(coeffs);
                % compute maximum of each sample
                dy = max(x, [], 1);
                % subtract approximated polynomial
                for j = 1:length(poolIds)
                    dy = dy - polyval(coeffs(j, :), x(j, :));
                end
                % compute maximum possible error

                L = interval(min(dy)-sum(tol), max(dy)+sum(tol));
                % add maximum possible error to output
                c_(outId) = c_(outId) + center(L);
                d(outId) = rad(L);
            end
        end
        % set output
        c = c_;
        G = G_;
        GI = GI_;
        [E, G] = removeRedundantExponents(E, G);

        % order reduction
        [c, G, GI, E, id, d] = nnHelper.reducePolyZono(c, G, GI, E, id, d, evParams.num_generators);

        % 4. add approximation error to polynomial zonotope

        c = c + 0.5 * d;

        % add error d
        Gd = diag(0.5*d);
        % Gd = diag(d);
        if evParams.add_approx_error_to_GI
            GI = [GI, Gd];
        else
            G = [G, Gd];
            E = blkdiag(E, eye(size(Gd, 2)));
            id = [id; 1 + (1:size(Gd, 2))' * id_];
            id_ = max(id);
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
