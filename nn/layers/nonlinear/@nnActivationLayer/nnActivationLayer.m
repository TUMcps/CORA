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

    function [nin, nout] = getNumNeurons(obj)
        nin = [];
        nout = [];
    end

    function outputSize = getOutputSize(obj, inputSize)
        outputSize = inputSize; % for most activation functions
    end
end

methods (Static)
    layer = instantiateFromString(activation)
end

methods (Abstract)
    [df_l,df_u] = getDerBounds(obj, l, u)
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
            if isempty(obj.l) || isempty(obj.u)
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

    function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
        % Obtain indices of active generators.
        genIds = obj.backprop.store.genIds;
        % Get size of generator matrix
        [n,q,batchSize] = size(G);
        % Compute radius of generators.
        r = reshape(sum(abs(G(:,genIds,:)),2),[n batchSize]);
        % Compute the bounds of the input.
        if options.nn.interval_center
            cl = reshape(c(:,1,:),[n batchSize]);
            cu = reshape(c(:,2,:),[n batchSize]);
            l = cl - r;
            u = cu + r;
        else
            l = c - r;
            u = c + r;
        end

        % Compute an image enclosure.
        [m,m_l,m_u,dl,dl_l,dl_u,du,du_l,du_u] = ...
            aux_imgEncBatch(obj,obj.f,obj.df,l,u,options,...
                @(m) obj.computeExtremePointsBatch(m,options));

        % Compute resulting generator matrix (without approx. errors).
        G = permute(m,[1 3 2]).*G;

        if options.nn.use_approx_error
            % Retrieve stored id-matrix and generator indices
            approxErrGenIds = obj.backprop.store.approxErrGenIds;
            % Retrieve number of approximation errors.
            if options.nn.interval_center
                dn = n;
            else
                dn = length(approxErrGenIds);
            end
            % Obtain number of generators after adding the approximation
            % errors.
            p = max(approxErrGenIds);
            % Identify for which dimensions to consider the approximation
            % error.
            % [~,dDims] = sort(1/2*(du - dl),1,'descend');
            % [~,dDims] = sort(rand([n batchSize],'like',c),1);
            dDims = repmat((1:n)',1,batchSize);
            dDimsIdx = reshape(sub2ind([n batchSize], ...
                dDims,repmat(1:batchSize,n,1)),size(dDims));
            notdDimsIdx = dDimsIdx(dn+1:end,:);
            dDimsIdx = dDimsIdx(1:dn,:);
            % % Set not considered approximation errors to 0.
            % dl(notdDimsIdx) = 0;
            % dl_l(notdDimsIdx) = 0;
            % dl_u(notdDimsIdx) = 0;
            % du(notdDimsIdx) = 0;
            % du_l(notdDimsIdx) = 0;
            % du_u(notdDimsIdx) = 0;
            % % Scale not-dropped dimensions.
            % if dn < n % && options.nn.train.noise > 0
            %     dropFac = n/dn;
            % else
            %     dropFac = 1;
            % end
            % % dropFac = 1;
            % m(dDimsIdx) = dropFac.*m(dDimsIdx);
            % dl(dDimsIdx) = dropFac.*dl(dDimsIdx);
            % du(dDimsIdx) = dropFac.*du(dDimsIdx);

            if ~options.nn.interval_center
                % Compute indices for approximation errors in the generator
                % matrix.
                GdIdx = reshape(sub2ind([n p batchSize], ...
                   reshape(dDims(1:dn,:),1,[]),...
                   repmat(approxErrGenIds,1,batchSize), ...
                   repelem(1:batchSize,1,dn)),[dn batchSize]);

                % Extend the generator matrix and add approximation errors.
               if q < p
                   % Append generators for the approximation errors
                   G = cat(2,G,zeros(n,length(approxErrGenIds),batchSize,'like',G));
               end
               % Add approximation errors to the generators.
               d = 1/2*(du - dl);
               G(GdIdx) = d(dDimsIdx);
            end
        end

        % Compute result.
        if options.nn.interval_center
            c = permute(cat(3,m.*cl - dl,m.*cu + du),[1 3 2]);
        else
            c = m.*c + 1/2*(du + dl);
        end

        % Store the gradients for backpropagation.
        if options.nn.train.backprop
            % Store the slope.
            obj.backprop.store.coeffs = m;
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
        % Obtain indices of active generators.
        genIds = obj.backprop.store.genIds;
        % Obtain stored slope.
        m = obj.backprop.store.coeffs;

        if options.nn.train.exact_backprop
            % Obtain the stored gradients.
            m_l = obj.backprop.store.m_l;
            m_u = obj.backprop.store.m_u;
            % Obtain the stored gradients.
            dl_l = obj.backprop.store.dl_l;
            dl_u = obj.backprop.store.dl_u;
            du_l = obj.backprop.store.du_l;
            du_u = obj.backprop.store.du_u;

            r_G = sign(G(:,genIds,:));
            m_G = permute(m_u - m_l,[1 3 2]).*r_G;

            [n,~,batchSize] = size(G);

            if options.nn.interval_center
                cl = reshape(c(:,1,:),[n batchSize]);
                cu = reshape(c(:,2,:),[n batchSize]);
                gl = reshape(gc(:,1,:),[n batchSize]);
                gu = reshape(gc(:,2,:),[n batchSize]);
                % Precompute outer product of gradients and inputs.
                hadProd_l = permute(gl.*cl,[1 3 2]);
                hadProd_u = permute(gu.*cu,[1 3 2]);
                hadProd = hadProd_l + hadProd_u + sum(gG(:,genIds,:).*G,2);

                % Backprop gradients.
                rgc_l = gl.*m + m_l.*reshape(hadProd,[n batchSize]);
                rgc_u = gu.*m + m_u.*reshape(hadProd,[n batchSize]);
                rgG = gG(:,genIds,:).*permute(m,[1 3 2]) + m_G.*hadProd; 

                if options.nn.use_approx_error
                    % Obtain indices of the approximation errors in the generator
                    % matrix.
                    dDimsIdx = obj.backprop.store.dDimsIdx;
                    notdDimsIdx = obj.backprop.store.notdDimsIdx;
                    % Compute gradients w.r.t. center and generators.
                    du_G = permute(du_u - du_l,[3 1 2]);
                    dl_G = permute(dl_u - dl_l,[3 1 2]);
    
                    % Backprop gradients.
                    rgc_l(dDimsIdx) = rgc_l(dDimsIdx) ...
                        + du_l.*gu(dDimsIdx) + dl_l.*gl(dDimsIdx);
                    % rgc_l(notdDimsIdx) = obj.df(c.inf(notdDimsIdx)).*gc.inf(notdDimsIdx);
                    rgc_u(dDimsIdx) = rgc_u(dDimsIdx) ...
                        + du_u.*gu(dDimsIdx) + dl_u.*gl(dDimsIdx);
                    % rgc_u(notdDimsIdx) = obj.df(c.sup(notdDimsIdx)).*gc.sup(notdDimsIdx);
    
                    rgG = permute(rgG,[2 1 3]);
                    r_G = permute(r_G,[2 1 3]);
                    rgG(:,dDimsIdx) = rgG(:,dDimsIdx) ... 
                        + (du_G(:,:).*gu(dDimsIdx(:))' ...
                            + dl_G(:,:).*gl(dDimsIdx(:))').*r_G(:,dDimsIdx);
                    % rgG(:,notdDimsIdx) = 0; % reshape(obj.df(c(notdDimsIdx)),1,[]).*gG(:,notdDimsIdx);
                    rgG = permute(rgG,[2 1 3]);     
                    % Assign results.
                    gc = permute(cat(3,rgc_l,rgc_u),[1 3 2]);
                    gG = rgG;
                else
                    % TODO
                end
            else
                m_c = (m_u + m_l);

                % Precompute outer product of gradients and inputs.
                hadProd = permute(gc.*c,[1 3 2]) + sum(gG(:,genIds,:).*G,2);
    
                % Backprop gradients.
                rgc = gc.*m + m_c.*reshape(hadProd,[n batchSize]);
                rgG = gG(:,genIds,:).*permute(m,[1 3 2]) + m_G.*hadProd;  

                if options.nn.use_approx_error
                    % Obtain indices of the approximation errors in the generator
                    % matrix.
                    GdIdx = obj.backprop.store.GdIdx;
                    dDimsIdx = obj.backprop.store.dDimsIdx;
                    notdDimsIdx = obj.backprop.store.notdDimsIdx;
                    % Compute gradients w.r.t. center and generators.
                    dc_c = 1/2*(du_u + du_l + dl_u + dl_l);
                    dc_G = 1/2*permute(du_u - du_l + dl_u - dl_l,[3 1 2]);
    
                    d_c = 1/2*(du_u + du_l - dl_u - dl_l);
                    d_G = 1/2*permute(du_u - du_l - dl_u + dl_l,[3 1 2]);
    
                    % Backprop gradients.
                    rgc(dDimsIdx) = rgc(dDimsIdx) ...
                        + dc_c.*gc(dDimsIdx) + d_c.*gG(GdIdx);
                    rgc(notdDimsIdx) = obj.df(c(notdDimsIdx)).*gc(notdDimsIdx);
    
                    rgG = permute(rgG,[2 1 3]);
                    r_G = permute(r_G,[2 1 3]);
                    rgG(:,dDimsIdx) = rgG(:,dDimsIdx) ... 
                        + dc_G(:,:).*r_G(:,dDimsIdx).*gc(dDimsIdx(:))' ...
                        + d_G(:,:).*r_G(:,dDimsIdx).*gG(GdIdx(:))';
                    % rgG(:,notdDimsIdx) = 0; % reshape(obj.df(c(notdDimsIdx)),1,[]).*gG(:,notdDimsIdx);
                    rgG = permute(rgG,[2 1 3]);     
                    % Assign results.
                    gc = rgc;
                    gG = rgG;
                else
                    % TODO
                end
            end
        else
            % Consider the approximation as fixed. Use the slope of the
            % approximation for backpropagation.
            if options.nn.interval_center
                gc = permute(m,[1 3 2]).*gc;
            else
                gc = gc.*m;
            end
            gG = gG(:,genIds,:).*permute(m,[1 3 2]);
        end
    end
end

% Auxiliary functions -----------------------------------------------------

methods

    function df_i = getDf(obj, i)
        df_i = nnHelper.lookupDf(obj,i);
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

    function [m,m_l,m_u,dl,dl_l,dl_u,du,du_l,du_u] = ...
            aux_imgEncBatch(obj,f,df,l,u,options,extremePoints)
        % Compute center and radius.
        c = 1/2*(u + l);
        r = 1/2*(u - l);
        % Compute slope.
        m = (f(u) - f(l))./(2*r);
        % Find indices where upper and lower bounds are equal.
        idxBoundsEq = abs(u - l) < eps('like',c); 
        % If lower and upper bound are too close, approximate the slope
        % at center; to prevent numerical issues.
        m(idxBoundsEq) = df(c(idxBoundsEq));

        % Compute gradient of the slope.
        if options.nn.train.backprop && ...
                options.nn.train.exact_backprop
            m_l = m./(2*r) - df(l)./(2*r);
            m_u = df(u)./(2*r) - m./(2*r);
            % Prevent numerical issues.
            ddf = obj.getDf(2);
            m_l(idxBoundsEq) = ddf(l(idxBoundsEq));
            m_u(idxBoundsEq) = ddf(u(idxBoundsEq));
        else
            % No gradient of the slope.
            m_l = 0;
            m_u = 0;
        end

        % Compute the approximation errors.
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
        else
            % No approximation errors. Use approximation errors for the
            % offset.
            dl = 1/2*(f(c) - m.*c);
            du = dl;
        end

        if options.nn.train.backprop && options.nn.train.exact_backprop
            if options.nn.use_approx_error
                if strcmp(options.nn.poly_method,'bounds')
                    % We only consider the lower bound. The approximation
                    % error at the lower and upper bound is equal.
                    x_l = cat(3,xs_m.*m_l,ones(size(l),'like',l));
                    x_u = cat(3,xs_m.*m_u,zeros(size(u),'like',u));
                else
                    x_l = cat(3,xs_m.*m_l,ones(size(l),'like',l), ...
                        zeros(size(l),'like',l));
                    x_u = cat(3,xs_m.*m_u,zeros(size(u),'like',u), ...
                        ones(size(u),'like',u));
                end
    
                % Compute gradient of the approximation errors.
                xl = xs(dlIdx);
                dfxlm = obj.df(xl) - m;
                dl_l = dfxlm.*x_l(dlIdx) - m_l.*xl;
                dl_u = dfxlm.*x_u(dlIdx) - m_u.*xl;
    
                xu = xs(duIdx);
                dfxum = obj.df(xu) - m;
                du_l = dfxum.*x_l(duIdx) - m_l.*xu;
                du_u = dfxum.*x_u(duIdx) - m_u.*xu;
            else
                dl_l = 1/2*(df(c) - m);
                dl_u = dl_l;
    
                du_l = dl_l;
                du_u = dl_l;
            end
        else
            % No gradients of the approximation errors.
            dl_l = 0;
            dl_u = 0;

            du_l = 0;
            du_u = 0;
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
