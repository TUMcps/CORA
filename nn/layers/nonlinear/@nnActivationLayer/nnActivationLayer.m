classdef (Abstract) nnActivationLayer < nnLayer
% nnActivationLayer - abstract class for non-linear layers
%
% Syntax:
%    obj = nnActivationLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Kochdumper, N., et al. (2022). Open-and closed-loop neural network
%        verification using polynomial zonotopes. 
%        arXiv preprint arXiv:2207.02715.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner, Lukas Koller
% Written:       28-March-2022
% Last update:   01-April-2022 (moved to class folder)
%                16-February-2023 (combined approx_type)
%                03-May-2023 (LK, backprop)
%                30-May-2023 (approx error/output bounds)
%                02-August-2023 (LK, zonotope batch-eval & -backprop)
%                19-August-2023 (zonotope batch-eval: memory optimizations for GPU training)
%                02-February-2024 (LK, better zonotope backpropagation)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = true     % whether the layer is refineable
end

properties
    % function handles

    f                       % function
    df                      % function derivative

    % adaptive refinement

    order = 1               % order of approximation polynomial
    refine_heu              % heuristic for refinement
    do_refinement = true    % whether the layer should be refined

    l = []                  % lower bound of last input
    u = []                  % upper bound of last input

    merged_neurons = []     % network reduction
end

methods
    % constructor
    function obj = nnActivationLayer(name)
        % call super class constructor
        obj@nnLayer(name)

        % init function handles
        obj.f = @(x) obj.evaluateNumeric(x, struct('backprop', false));
        obj.df = obj.getDf(1);
    end
end

% evaluate (element-wise) -------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
       S = S.*permute(obj.df(x),[3 1 2]);
    end

    % interval
    function bounds = evaluateInterval(obj, bounds, options)
        if options.nn.reuse_bounds
            % save bounds
            if isempty(obj.l) || isempty(obj.u) || ~all(size(bounds) == size(obj.l))
                obj.l = bounds.inf;
                obj.u = bounds.sup;

                % set bounds
            elseif representsa_(bounds,'emptySet',eps) && ...
                    any(isnan(obj.l)) && any(isnan(obj.u))
                bounds = interval(obj.l, obj.u);
            end

            obj.l = max(obj.l, bounds.inf);
            obj.u = min(obj.u, bounds.sup);
        end

        % propagate through layer
        bounds = evaluateInterval@nnLayer(obj, bounds, options);
    end

    % zonotope/polyZonotope
    [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
    [c, G, GI, d] = evaluatePolyZonotopeNeuron(obj, c, G, GI, E, Es, order, ind, ind_, options)

    function [rc, rG] = evaluateZonotopeBatch(obj, c, G, options)
        % Compute image enclosure
        [rc,rG,coeffs] = aux_imgEncBatch(obj,obj.f,obj.df,c,G,options, ...
            @(m) obj.computeExtremePointsBatch(m,options));

        % store inputs and coeffs for backpropagation
        if options.nn.train.backprop
            % Store coefficients
            obj.backprop.store.coeffs = coeffs;

            % Store the slope.
            if options.nn.train.exact_backprop
                % Store gradient for the backprop through an image
                % enclosure.
                obj.backprop.store.m_l = m_l;
                obj.backprop.store.m_u = m_u;

                if options.nn.use_approx_error
                    if ~options.nn.interval_center
                        obj.backprop.store.GdIdx = GdIdx;
                    end
                    obj.backprop.store.dDimsIdx = dDimsIdx;
                    obj.backprop.store.notdDimsIdx = notdDimsIdx;

                    obj.backprop.store.dl_l = dl_l(dDimsIdx);
                    obj.backprop.store.dl_u = dl_u(dDimsIdx);
                    obj.backprop.store.du_l = du_l(dDimsIdx);
                    obj.backprop.store.du_u = du_u(dDimsIdx);
                end
            end
        end
    end
    
    % taylm
    r = evaluateTaylm(obj, input, options)
    function r = evaluateTaylmNeuron(obj, input, order, options)
        % enclose the ReLU activation function with a Taylor model by
        % fitting a quadratic function

        % compute lower and upper bound
        int = interval(input);
        l = infimum(int);
        u = supremum(int);

        % compute approx poly + error
        [coeffs, d] = computeApproxPoly(obj, l, u, order, options.nn.poly_method);

        % evaluate
        r = coeffs(end) + interval(-d, d);
        for i=1:length(coeffs)-1
            r = r + coeffs(end-i) * input^i;
        end        
    end
    
    % conZonotope
    [c, G, C, d, l, u] = evaluateConZonotope(obj, c, G, C, d, l, u, options)
    function [c, G, C, d, l, u] = evaluateConZonotopeNeuron(obj, c, G, C, d, l, u, j, options)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'conZonotope'))
    end

    % backprop ------------------------------------------------------------

    function grad_in = backpropNumeric(obj, input, grad_out, options)
        % backpropagte gradient
        grad_in = obj.df(input) .* grad_out;
    end

    function [gl, gu] = backpropIntervalBatch(obj, l, u, gl, gu, options)
        gl = obj.df(l) .* gl;
        gu = obj.df(u) .* gu;
    end

    function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options)
        % obtain stored coefficients of image enclosure from forward prop.
        coeffs = obj.backprop.store.coeffs;
        % obtain slope of the approximation
        m = reshape(coeffs(:,1,:),size(coeffs,[1 3]));
        % obtain indices of active generators
        genIds = obj.backprop.store.genIds;
        if options.nn.train.exact_backprop
            % Obtain indices of the approximation errors in the generator
            % matrix.
            GdIdx = obj.backprop.store.GdIdx;
            dDimsIdx = obj.backprop.store.dDimsIdx;

            m_c = obj.backprop.store.m_c;
            m_G = obj.backprop.store.m_G;

            dc_c = obj.backprop.store.dc_c;
            dc_G = obj.backprop.store.dc_G;

            d_c = obj.backprop.store.d_c;
            d_G = obj.backprop.store.d_G;

            % Precompute outer product of gradients and inputs.
            hadProdc = permute(gc.*c,[1 3 2]);
            hadProdG = gG(:,genIds,:).*G;

            % Backprop gradients.
            rgc = gc.*m + m_c.*reshape(hadProdc + sum(hadProdG,2),size(c));
            rgc(dDimsIdx) = rgc(dDimsIdx) + dc_c.*gc(dDimsIdx);
            rgc(dDimsIdx) = rgc(dDimsIdx) + d_c.*gG(GdIdx);
            % Assign results.
            gc = rgc;

            rgG = gG(:,genIds,:).*permute(m,[1 3 2]) + m_G.*(hadProdc + hadProdG);
            rgG = permute(rgG,[2 1 3]);
            rgG(:,dDimsIdx) = rgG(:,dDimsIdx) + dc_G.*reshape(gc(dDimsIdx),1,[]);
            rgG(:,dDimsIdx) = rgG(:,dDimsIdx) + d_G.*reshape(gG(GdIdx),1,[]);
            rgG = permute(rgG,[2 1 3]);          
            % Assign results.
            gG = rgG;

        else
            % Consider the approximation as fixed. Use the slope of the
            % approximation for backpropagation
            gc = gc.*m;
            gG = gG(:,genIds,:).*permute(m,[1 3 2]);
        end

        % Clear backprop storage.
        clear 'obj.backprop.store';
    end
end

% Auxiliary functions -----------------------------------------------------

methods

    function df_i = getDf(obj, i)
        df_i = nnHelper.lookupDf(obj,i);
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = inputSize; % for most activation functions
    end

    % approximation polynomial + error

    [coeffs, d] = computeApproxPoly(obj, l, u, varargin)

    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        % bound approximation error according to [1, Sec. 3.2]

        % compute the difference between activation function and quad. fit
        [df_l,df_u] = obj.getDerBounds(l, u);
        [diffl,diffu] = nnHelper.minMaxDiffOrder(coeffs, l, u, obj.f, df_l,df_u);
        
        % change polynomial s.t. lower and upper error are equal
        diffc = (diffl+diffu)/2;
        coeffs(end) = coeffs(end) + diffc;
        d = diffu-diffc; % error is radius then.
    end

    % find regions with approximating polynomials
    coeffs = findRegionPolys(obj,tol,order,l_max,u_max,pStart,dStart,pEnd,dEnd)

    function fieldStruct = getFieldStruct(obj)
        fieldStruct = struct;
        if ~isempty(obj.merged_neurons)
            fieldStruct.merged_neurons = obj.merged_neurons;
        end
    end
end

methods (Static)
    layer = instantiateFromString(activation)
end

methods (Abstract)
    [df_l,df_u] = getDerBounds(obj, l, u)
end

methods(Access=protected)
    function [coeffs, d] = computeApproxPolyCustom(obj, l, u, order, poly_method)
        % implement custom polynomial computation in subclass
        coeffs = []; d = [];
    end

    function [xs,xs_m] = computeExtremePointsBatch(obj, m, options)
        throw(CORAerror('CORA:nnLayerNotSupported', obj, 'computeExtremePointsBatch'))
    end

    % function xs = computeExtremePointsDfBatch(obj, m)
    %     throw(CORAerror('CORA:nnLayerNotSupported', obj, 'computeExtremePointsDfBatch'))
    % end

    function [rc,rG,coeffs,d] = aux_imgEncBatch(obj,f,df,c,G,options,extremePoints)
        % obtain indices of active generators
        genIds = obj.backprop.store.genIds;
        % compute bounds
        r = reshape(sum(abs(G(:,genIds,:)),2),size(c));
        % r = max(eps('like',c),r); % prevent numerical instabilities
        l = c - r;
        u = c + r;
        % compute slope of approximation
        if strcmp(options.nn.poly_method,'bounds')
            m = (f(u) - f(l))./(2*r);
            % indices where upper and lower bound are equal
            idxBoundsEq = abs(u - l) < eps('like',c); 
            % If lower and upper bound are too close, approximate the slope
            % at center.
            m(idxBoundsEq) = df(c(idxBoundsEq));
            if options.nn.train.backprop
                obj.backprop.store.idxBoundsEq = idxBoundsEq;
            end
        elseif strcmp(options.nn.poly_method,'center')
            m = cast(df(c),'like',c);
        elseif strcmp(options.nn.poly_method,'singh')
            m = min(cast(df(l),'like',l),cast(df(u),'like',u));
        else
            throw(CORAerror('CORA:wrongFieldValue', ...
                sprintf("Unsported 'options.nn.poly_method': %s",...
                    options.nn.poly_method)));
        end

        % evaluate image enclosure
        rc = m.*c;
        % rG(:,genIds,:) = permute(m,[1 3 2]).*G(:,genIds,:);
        rG = permute(m,[1 3 2]).*G;

        if options.nn.use_approx_error
           % Compute extreme points.
           [xs,xs_m] = extremePoints(m);
           % Determine number of extreme points.
           s = size(xs,3);
           % Add interval bounds.
           if strcmp(options.nn.poly_method,'bounds')
               % the approximation error at l and u are equal, thus we only
               % consider the upper bound u.
               xs = cat(3,xs,l);
           else
               xs = cat(3,xs,l,u);
           end
           ys = f(xs);
           % Compute approximation error at candidates.
           ds = ys - m.*xs;
           % We only consider candidate extreme points within boundaries.
           notInBoundsIdx = (xs < l | xs > u);
           ds(notInBoundsIdx) = inf;
           [dl,dlIdx] = min(ds,[],3,'linear');
           ds(notInBoundsIdx) = -inf;
           [du,duIdx] = max(ds,[],3,'linear');

           % Retrieve stored id-matrix and generator indices
           approxErrGenIds = obj.backprop.store.approxErrGenIds;
           % Retrieve number of approximation errors.
           dn = length(approxErrGenIds);
           % Get size of generator matrix
           [n,q,batchSize] = size(rG);
           p = max(approxErrGenIds);
           if q < p
               % Append generators for the approximation errors
               rG = cat(2,rG,zeros(n,length(approxErrGenIds),batchSize,'like',rG));
           end
           % Obtain the dn largest approximation errors.
           % [~,dDims] = sort(1/2*(du - dl),1,'descend');
           dDims = repmat((1:n)',1,batchSize);
           dDimsIdx = reshape(sub2ind([n batchSize], ...
               dDims,repmat(1:batchSize,n,1)),size(dDims));
           notdDimsIdx = dDimsIdx(dn+1:end,:);
           % set not considered approx. error to 0
           dl(notdDimsIdx) = f(c(notdDimsIdx)) - m(notdDimsIdx).*c(notdDimsIdx);
           du(notdDimsIdx) = dl(notdDimsIdx);
           % shift y-intercept by center of approximation errors
           t = 1/2*(du + dl);
           d = 1/2*(du - dl);

           % Compute indices for approximation errors in the generator
           % matrix.
           GdIdx = reshape(sub2ind([n p batchSize], ...
               repmat(1:n,1,batchSize), ...reshape(dDims(1:dn,:),1,[]),...
               repmat(approxErrGenIds,1,batchSize), ...
               repelem(1:batchSize,1,n)),[dn batchSize]);
           % Store indices of approximation error in generator matrix.
           obj.backprop.store.GdIdx = GdIdx;
           dDimsIdx = dDimsIdx(1:dn,:);
           obj.backprop.store.dDimsIdx = dDimsIdx;
           % Add approximation errors to the generators.
           rG(GdIdx) = d; % (dDimsIdx);

           % compute gradients
           if options.nn.train.backprop && options.nn.train.exact_backprop
               % derivative of radius wrt. generators
               r_G = sign(G);
               % Compute gradient of the slope.
               if strcmp(options.nn.poly_method,'bounds')
                   m_c = (obj.df(u) - obj.df(l))./(2*r);
                   m_G = r_G.*permute((obj.df(u) + obj.df(l) - 2*m)./(2*r),[1 3 2]);
                   % prevent numerical issues
                   ddf = obj.getDf(2);
                   m_c(idxBoundsEq) = obj.df(c(idxBoundsEq));
                   m_G = permute(m_G,[2 1 3]);
                   r_G = permute(r_G,[2 1 3]);
                   m_G(:,idxBoundsEq) = r_G(:,idxBoundsEq).*ddf(c(idxBoundsEq))';
                   m_G = permute(m_G,[2 1 3]);
                   r_G = permute(r_G,[2 1 3]);
               elseif strcmp(options.nn.poly_method,'center')
                   ddf = obj.getDf(2);
                   m_c = ddf(inc).*ones(size(inc),'like',c); 
                   m_G = zeros(size(inG),'like',inG);
               elseif strcmp(options.nn.poly_method,'singh')
                   lu = cat(3,l,u);
                   [~,mIdx] = min(obj.df(lu),[],3,'linear');
                   ddf = obj.getDf(2);
                   m_c = ddf(lu(mIdx)).*ones(size(c),'like',c);
                   [~,~,mIdx] = ind2sub(size(lu),mIdx);
                   m_G = permute(-1*(mIdx == 1).*m_c,[1 3 2]).*drdG;
               else
                   throw(CORAerror('CORA:wrongFieldValue', ...
                       sprintf("Unsported 'options.nn.poly_method': %s",...
                           options.nn.poly_method)));
               end

               % Add gradients for interval bounds.
               if strcmp(options.nn.poly_method,'bounds')
                   % We only consider the lower bound. The approximation
                   % error at the lower and upper bound is equal.
                   x_c = cat(3,xs_m.*m_c,ones(size(l),'like',l));
                   x_G = cat(4,permute(xs_m,[1 2 4 3]).*permute(m_G,[1 3 2]), ...
                       -permute(r_G,[1 3 2]));
               else
                   x_c = cat(3,xs_m.*mc,ones(size(l),'like',l), ...
                       ones(size(u),'like',u));
                   x_G = cat(4,permute(xs_m,[1 2 4 3]).*permute(m_G,[1 3 2]), ...
                       -permute(r_G,[1 3 2]),permute(r_G,[1 3 2]));
               end
               x_G = permute(x_G,[3 1 2 4]);

               % Compute gradient of the approximation errors.
               xl = xs(dlIdx);
               dfxlm = obj.df(xl) - m;
               dl_c = dfxlm.*x_c(dlIdx) - m_c.*xl;
               xl_G = reshape(x_G(:,dlIdx),[q n batchSize]);
               dl_G = permute(dfxlm,[3 1 2]).*xl_G - permute(m_G,[2 1 3]).*permute(xl,[3 1 2]);

               xu = xs(duIdx);
               dfxum = obj.df(xu) - m;
               du_c = dfxum.*x_c(duIdx) - m_c.*xu;
               xu_G = reshape(x_G(:,duIdx),[q n batchSize]);
               du_G = permute(dfxum,[3 1 2]).*xu_G - permute(m_G,[2 1 3]).*permute(xu,[3 1 2]);

               % Compute components for the backpropagation.
               obj.backprop.store.m_c = m_c;
               obj.backprop.store.m_G = m_G;

               obj.backprop.store.dc_c = 1/2*(du_c(dDimsIdx) + dl_c(dDimsIdx));
               obj.backprop.store.dc_G = ...
                   1/2*(du_G(:,dDimsIdx) + dl_G(:,dDimsIdx));

               obj.backprop.store.d_c = 1/2*(du_c(dDimsIdx) - dl_c(dDimsIdx));
               obj.backprop.store.d_G = ...
                   1/2*(du_G(:,dDimsIdx) - dl_G(:,dDimsIdx));

           end
        else
            % compute y-intercept
            t = f(c) - m.*c;
            % approximation errors are 0
            d = 0;
        end
        % return coefficients
        coeffs = permute(cat(3,m,t),[1 3 2]);
        % Add y-intercept.
        rc = rc + t;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
