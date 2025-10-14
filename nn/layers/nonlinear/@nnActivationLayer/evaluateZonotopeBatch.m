function [c, G] = evaluateZonotopeBatch(obj, c, G, options)
% evaluateZonotopeBatch - evaluate nnActivationLayer for a batch of 
%   zonotopes
%
% Syntax:
%    [c, G] = layeri.evaluateZonotopeBatch(c, G, options);
%
% Inputs:
%    c, G - batch of zonotope; [n,q+1,b] = size([c G]),
%       where n is the number of dims, q the number of generators, and b the batch size
%    options - parameter for neural network evaluation
%
% Outputs:
%    c, G - batch of output sets
%
% References:  
%    [1] Koller et al. Set-based Training for Neural Network Verification.
%       (TMLR). 2025
%    [2] Singh et al. Fast and Effective Robustness Certification.
%       (NeurIPS). 2018
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/evaluateZonotopeBatch

% Authors:       Lukas Koller
% Written:       12-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Obtain indices of active generators.
genIds = obj.genIds;
% Get size of generator matrix.
[n,q,bSz] = size(G);

% Compute radius of input sets.
r = reshape(sum(abs(G(:,genIds,:)),2),[n bSz]);
% Compute the bounds of the center.
if options.nn.interval_center
    % The center store the approximation errors.
    cl = reshape(c(:,1,:),[n bSz]);
    cu = reshape(c(:,2,:),[n bSz]);
else
    % The center is just a single point.
    cl = c;
    cu = c;
end
% Compute the bounds of the input.
l = cl - r;
u = cu + r;

if isfield(obj.backprop.store,'l') ...
    && all(size(l) == size(obj.backprop.store.l)) ...
    && isfield(obj.backprop.store,'u') ...
    && all(size(u) == size(obj.backprop.store.u))
    % Obtain stored bounds.
    l = max(l,min(u,obj.backprop.store.l));
    u = min(u,max(l,obj.backprop.store.u));
end

% Compute an image enclosure.
[m,m_l,m_u,dl,dl_l,dl_u,dl_m,du,du_l,du_u,du_m] = ...
    aux_imgEncBatch(obj,obj.f,obj.df,l,u,options,...
        @(m) obj.computeExtremePointsBatch(m,options));

% Compute resulting generator matrix (without approx. errors).
G = permute(m,[1 3 2]).*G;

% Check if approximation error should be considered; some set-training
% schemes ignore the approximation errors.
if options.nn.use_approx_error
    % Retrieve the indices of the generators for the approximation errors.
    % The indices are stored before the evaluation with a call of
    % neuralNetwork/prepareForZonoBatchEval.
    approxErrGenIds = obj.approxErrGenIds;
    % Obtain the number of considered approximation errors.
    dn = length(approxErrGenIds);
    % Obtain number of generators after adding the approximation
    % errors.
    p = max(approxErrGenIds);
    % Identify for which dimensions to consider the approximation
    % error.
    if n == dn
        % There is enough space for all approximation errors.
        options.nn.approx_error_order = 'sequential';
    end
    switch options.nn.approx_error_order
        case 'length'
            % Use the dimensions with the largest intervals.
            [~,dDims] = sort(1/2*(du - dl),1,'descend');
        case 'sensitivity*length'
            % Use the most sensitive dimensions with the largest intervals.

            % The sensitivity has to be stored beforehand with
            % neuralNetwork/calcSensitivity and storeSensitivity=true.
            % Check if the sensitivity was stored.
            if ~isprop(obj,'sensitivity')
                throw(CORAerror('CORA:notSupported', ...
                    ['The option options.nn.approx_error_order=' ...
                    '"sensitivity*length" can only be used if the ' ...
                    'sensitivity was stored beforehand with a call of ' ...
                    'neuralNetwork/calcSensitivity and ' ...
                    'storeSensitivity=true']));
            end
            % Check if the stored sensitivity has the correct batch size.
            % If not, there was neuron splitting involved by
            % neuralNetwork/verify. We match the batch size.
            if size(obj.sensitivity,3) ~= bSz
                newSplits = bSz/size(obj.sensitivity,3);
                obj.sensitivity = repelem(obj.sensitivity,1,1,newSplits);
            end
            % Obtains sensitivity matrix and sum across outputs.
            S = reshape(sum(abs(obj.sensitivity),1),[n bSz]);
            % Find dimensions with the hightest heuristic.
            [~,dDims] = sort(1/2*(du - dl).*S,1,'descend');
        case 'random'
            % Select random dimensions.
            [~,dDims] = sort(rand([n bSz],'like',c),1);
        case 'sequential'
            % Enumerate the dimensions sequentially.
            dDims = repmat((1:n)',1,bSz);
        otherwise
    end
    % Compute the indices into for the considered approximation errors.
    dDimsIdx = reshape(sub2ind([n bSz],dDims, ...
        repmat(1:bSz,n,1)),size(dDims));
    % Extract the indices for approximation errors.
    dDimsIdx = dDimsIdx(1:dn,:);
    % The remaining indices are for the approximation errors that are not
    % considered.
    notdDimsIdx = dDimsIdx(dn+1:end,:);
    % If there are considered approximation errors we add them into the
    % generator matrix.
    if dn > 0
        % We might need to extend the generator matrix to add the 
        % approximation errors.
        if q < p
            % Append sufficiently many generators for the approximation 
            % errors.
            G = [G zeros([n p-q bSz],'like',G)];
        end
        % Compute indices for approximation errors into the generator
        % matrix.
        GdIdx = reshape(sub2ind(size(G), ...
           reshape(dDims(1:dn,:),1,[]),...
           repmat(approxErrGenIds,1,bSz), ...
           repelem(1:bSz,1,dn)),[dn bSz]);

        % Compute the considered approximation errors.
        d = 1/2*(du(dDimsIdx) - dl(dDimsIdx));
        % Use the computed indices to write the approximation errors into 
        % the correct generators.
        G(GdIdx) = d;
    end
end

% Compute the resulting center.
if options.nn.interval_center
    % Compute offset for the center.
    offset = 1/2*(du + dl);
    % The offset is only applied to dimensions for which the approximation
    % error is stored in the generator matrix.
    offset(notdDimsIdx) = 0;
    % Set the considered approximation errors to 0; only the remaining
    % approximation error are added to the center interval.
    dcl = dl;
    dcl(dDimsIdx) = 0;
    dcu = du;
    dcu(dDimsIdx) = 0;
    % The center store the remaining approximation errors.
    c = permute(cat(3,m.*cl + offset - dcl,m.*cu + offset + dcu),[1 3 2]);
else
    % We apply the computed slope to the center and add an offset for the
    % approximation errors.
    c = m.*c + 1/2*(du + dl);
end

% Store the gradients for backpropagation.
if options.nn.train.backprop
    % Store the slope.
    obj.backprop.store.coeffs = m;
    % Store the approximation erros.
    obj.backprop.store.dl = dl;
    obj.backprop.store.du = du;
    % The flag exact_backprop toggles the exact backpropagation through the
    % image enclosure. For that we have to store additional gradients.
    if options.nn.train.exact_backprop
        % Store the gradient of the slope w.r.t. the input bounds l and u.
        obj.backprop.store.m_l = m_l;
        obj.backprop.store.m_u = m_u;
        % Store the gradients of the approximation errors w.r.t. the input
        % bounds l and u; additionally, store index information.
        if options.nn.use_approx_error
            % obj.backprop.store.GdIdx = GdIdx;
            obj.backprop.store.dDimsIdx = dDimsIdx;
            obj.backprop.store.notdDimsIdx = notdDimsIdx;
            obj.backprop.store.dl_l = dl_l;
            obj.backprop.store.dl_u = dl_u;
            obj.backprop.store.dl_m = dl_m;
            obj.backprop.store.du_l = du_l;
            obj.backprop.store.du_u = du_u;
            obj.backprop.store.du_m = du_m;
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [m,m_l,m_u,dl,dl_l,dl_u,dl_m,du,du_l,du_u,du_m] = ...
        aux_imgEncBatch(obj,f,df,l,u,options,extremePoints)
    % Compute center and radius of the input set.
    c = 1/2*(u + l);
    r = 1/2*(u - l);
    switch options.nn.poly_method
        case 'bounds'
            % Compute slope.
            m = (f(u) - f(l))./(2*r);

            if options.nn.verif_slope_optim_step_size > 0 ...
                    && isfield(obj.backprop.store,'dm') ...
                    && all(size(m) == size(obj.backprop.store.dm))
                % Gradient-based optimization of the slope. Cannot be 
                % negative; most activation function are monotonic.
                m = m + obj.backprop.store.dm; 
            end

            % Compute gradient of the slope.
            if options.nn.train.backprop && options.nn.train.exact_backprop
                m_l = (-df(l) + m)./(2*r);
                m_u = (df(u) - m)./(2*r);
            end
        % TODO: fix the other cases.
        case 'singh' % Not supported for ReLU
            % See [2].
            [lambda,idx] = min([df(l) df(u)],[],2);
            m = (f(u) - f(l) - lambda.*(2*r)).*(2*r);
            % Compute gradient of the slope.
            if options.nn.train.backprop && options.nn.train.exact_backprop
                ddf = obj.getDf(2);
                dlambda_l = [ddf(l).*(2*r) - lambda, -df(u)];
                m_l = df(l) + dlambda_l(idx);
                dlambda_u = [-df(l), ddf(u).*(2*r) - lambda];
                m_u = df(u) - dlambda_u(idx);
            end
        case 'center'
            % Compute slope.
            m = df(c);
            % Compute gradient of the slope.
            if options.nn.train.backprop && options.nn.train.exact_backprop
                ddf = obj.getDf(2);
                m_l = 1/2*ddf(c);
                m_u = 1/2*ddf(c);
            end
        otherwise
            throw(CORAerror('CORA:wrongFieldValue',...
                'options.nn.poly_method',{'bounds','singh','center'}));
    end

    % Find indices where upper and lower bounds are equal.
    idxBoundsEq = withinTol(u,l,eps('like',c)); 
    % If lower and upper bound are too close, approximate the slope
    % at center; to prevent numerical issues.
    m(idxBoundsEq) = df(c(idxBoundsEq));
    % Prevent numerical issues.
    if options.nn.train.backprop && ...
            options.nn.train.exact_backprop
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

        % For 'bounds' without slope optimization the approximation error at
        % the lower and the upper bound are equal.
        lowerAndUpperErrorEqual = strcmp(options.nn.poly_method,'bounds') ...
            && options.nn.verif_slope_optim_step_size == 0;

        % Add interval bounds.
        if lowerAndUpperErrorEqual
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
            if lowerAndUpperErrorEqual
                % We only consider the lower bound. The approximation
                % error at the lower and upper bound is equal.
                x_l = cat(3,xs_m.*m_l,ones(size(l),'like',l));
                x_u = cat(3,xs_m.*m_u,zeros(size(u),'like',u));

                x_m = cat(3,xs_m,zeros(size(l),'like',l));
            else
                x_l = cat(3,xs_m.*m_l,ones(size(l),'like',l), ...
                    zeros(size(u),'like',u));
                x_u = cat(3,xs_m.*m_u,zeros(size(l),'like',l), ...
                    ones(size(u),'like',u));

                x_m = cat(3,xs_m,zeros(size(l),'like',l), ...
                    ones(size(u),'like',u));
            end

            % Compute gradient of the approximation errors.
            xl = xs(dlIdx);
            dfxlm = obj.df(xl) - m;
            dl_l = dfxlm.*x_l(dlIdx) - m_l.*xl; % grad of dl w.r.t. l
            dl_u = dfxlm.*x_u(dlIdx) - m_u.*xl; % grad of dl w.r.t. u

            xu = xs(duIdx);
            dfxum = obj.df(xu) - m;
            du_l = dfxum.*x_l(duIdx) - m_l.*xu; % grad of du w.r.t. l
            du_u = dfxum.*x_u(duIdx) - m_u.*xu; % grad of du w.r.t. u

            % The gradient of the approximation errors w.r.t. the slope.
            dl_m = obj.df(xl).*x_m(dlIdx) - xl;
            du_m = obj.df(xu).*x_m(duIdx) - xu;
        else
            % The center is only shifted by the 1/2*(du + dl); thus, the
            % gradient are all the same.
            dl_l = 1/2*(df(c) - m);
            dl_u = dl_l;

            du_l = dl_l;
            du_u = dl_l;

            % The gradient of the approximation errors w.r.t. the slope is 0.
            dl_m = 0;
            du_m = 0;
        end
    else
        % No gradients of the approximation errors.
        dl_l = 0;
        dl_u = 0;
        du_l = 0;
        du_u = 0;
        dl_m = 0;
        du_m = 0;
    end
end

% ------------------------------ END OF CODE ------------------------------
