function [l,u,bl,bu] = approximateBoundsWithGPU(cZs,numUnionConst,options)
% approximateBoundsWithGPU - Approximate the bounds of a batched
%   constrained zonotope. Mainly, used for the input refinement procedure 
%   in neuuralNetwork/verify. This function is tuned towards speed rather
%   than accuracy; it significantly benefits from GPU computations. 
%   This is an efficent implementation of [1, Prop 3]. 
%   WARNING: The input zonotopes are constraint zonotopes with inequality
%       constraints; not equality- Please transform accordingly!
%
% Syntax:
%    cZs.c = <center [n x batch size]>;
%    cZs.G = <generator matrix [n x q x batch size]>; % q: number of generators.
%    cZs.A = <constraint matrix (on factor space) [p x q x batch size]>; % p: number of constraints.
%    cZs.b = <constraint offset [p x batch size]>;
%    cZs.dr = <interval error [n x batch size]>; % optional interval error
%       on the zonotope.
%    numUnionConst = 1; % Compute the intersection of all constraints.
%    options.nn.polytope_bound_approx_max_iter = 4; % Maximum number of iterations
%    options.nn.exact_conzonotope_bounds = false; % Approximate the bounds.
%    options.nn.batch_union_conzonotope_bounds = true; % Move union constraint to the batch.
%    [l,u,bl,bu] = approximateBoundsWithGPU(cZs,numUnionConst,options)
%
% Inputs:
%    cZs - (struct) batch constraint zonotope, please see the syntax for 
%       an example. For GPU ultilization, move everything to GPU before,
%       i.e., x = cast(x,'single gpuArray');
%    numUnionConst - number of unions constraints; the first 
%       #numUnionConst constraints of cZs are unified (needed for 
%       safeSet specifications).
%    options.nn.conzonotope_bounding_method - bounding method, i.e., 
%       {'fourier-motzkin','dual-iter','exact'}
%    options.nn.batch_union_conzonotope_bounds - (bool) batch union 
%       constraints
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Koller, L. "Out of the Shadows: Exploring a Latent Space for 
%       Neural Network Verification". ICLR (2026)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/verify

% Authors:       Lukas Koller
% Written:       19-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Input arguments represent a constraint zonotope with inequality
% constraints.
% numUnionConst: number of unions constraints; the first #numUnionConst
% constraints of cZs are unified (needed for safeSet specifications).
% options.nn.conzonotope_bounding_method: use linear programs to compute the bounds.
% options.nn.batch_union_conzonotope_bounds: batch union constraints

% Extract parameters of the constraint zonotope.
c = cZs.c;
G = cZs.G;
dr = cZs.dr;
A = cZs.A;
b = cZs.b;

% Obtain number of dimensions, generators, and batch size.
[n,q,bSz] = size(G);

if isempty(A)
    % There are no constraints. Just compute the bounds of the
    % zonotope.
    r = dr + reshape(sum(abs(G),2),[n bSz]);
    l = c - r;
    u = c + r;
    % The bounds of the hypercube are just -1 and 1;
    bl = -ones([q bSz],'like',G);
    bu = ones([q bSz],'like',G);
    return;
end

% Pad the generator matrix with 0s if there are more dimensions in the 
% hypercube.
G = cat(2,G,zeros([n size(A,2)-q bSz],'like',G));
% Update the number of generators.
[~,q,~] = size(G);

% Specify indices of intersection constraints.
intConIdx = (numUnionConst+1):size(A,1);

% The safe set is the union of all constraints. Thus, we have to create a 
% new set for each constraint. Either compute the bounds in batch-wise
% fashion (for speed) or for-loop (to save memory space).
if options.nn.batch_union_conzonotope_bounds
    % Move union constraints into the batch.
    Au = reshape(permute(A(1:numUnionConst,:,:),[4 2 3 1]),...
        [1 q bSz*numUnionConst]);
    bu = reshape(permute(b(1:numUnionConst,:),[3 2 1]),...
        [1 bSz*numUnionConst]);
    % Replicate intersection constraints.
    Ai = repmat(A(intConIdx,:,:),1,1,numUnionConst);
    bi = repmat(b(intConIdx,:),1,numUnionConst);
    % Append intersection constraints.
    A = cat(1,Au,Ai);
    b = cat(1,bu,bi);

    % Replicate the center and generator matrix.
    c = repmat(c,1,numUnionConst);
    if numel(dr) > 1
        dr = repmat(dr,1,numUnionConst);
    end
    G = repmat(G,1,1,numUnionConst);

    % Approximate the bounds of the hypercube (bounded polytope).
    [l,u,bl,bu] = aux_boundConZonotope(c,G,dr,A,b,options);

    if numUnionConst > 1
        % Unify sets if a safe set is specified.
        l = min(reshape(l,[n bSz numUnionConst]),[],3);
        u = max(reshape(u,[n bSz numUnionConst]),[],3);
        bl = min(reshape(bl,[q bSz numUnionConst]),[],3);
        bu = max(reshape(bu,[q bSz numUnionConst]),[],3);
    end
else
    % Initialize the bounds.
    l = inf([n bSz],'like',G);
    u = -inf([n bSz],'like',G);
    bl = -ones([q bSz],'like',G);
    bu = ones([q bSz],'like',G);

    % Loop over the union constraints.
    for k=1:numUnionConst
        % Use the k-th union constraint and all intersection
        % constraints.
        Ak = A([k intConIdx],:,:);
        bk = b([k intConIdx],:);
        % Approximate the bounds of the hypercube.
        [lk,uk,blk,buk] = aux_boundConZonotope(c,G,dr,Ak,bk,options);
        % Unify constraints.
        l = min(l,lk);
        u = max(u,uk);
        bl = min(bl,blk);
        bu = max(bu,buk);
    end
end

% Identify empty sets.
isEmpty = any(bl > bu,1) | any(isnan(l) | isnan(u),1);
l(:,isEmpty) = NaN;
u(:,isEmpty) = NaN;
bl(:,isEmpty) = 0;
bu(:,isEmpty) = 0;

end


% Auxiliary functions -----------------------------------------------------

function [l,u,bl,bu] = aux_boundConZonotope(c,G,dr,A,b,options)
    % Compute the bounds [l,u] of a constrained zonotope.

    % Obtain number of dimensions, generators, and batch size.
    [n,q,bSz] = size(G);

    switch options.nn.conzonotope_bounding_method
        case 'fourier-motzkin'
            % Approximate the bounds of the costrained hypercube.
            [bl,bu] = aux_fourierMotzkinApproximation(A,b,options);

            % Map bounds of the factors to bounds of the constraint 
            % zonotope. We use interval arithmetic for that.
            bc = 1/2*permute(bu + bl,[1 3 2]);
            br = 1/2*permute(bu - bl,[1 3 2]);
            % Map bounds of the factors to bounds of the constraint 
            % zonotope.
            c = c + reshape(pagemtimes(G,bc),[n bSz]);
            r = dr + reshape(pagemtimes(abs(G),br),[n bSz]);
            l = c - r;
            u = c + r;
        case 'dual-iter'
            % Approximate the bounds of the costrained hypercube for a
            % simple feasibility check.
            [bl,bu] = aux_fourierMotzkinApproximation(A,b,options);

            % Map bounds of the factors to bounds of the constraint 
            % zonotope. We use interval arithmetic for that.
            bc = 1/2*permute(bu + bl,[1 3 2]);
            br = 1/2*permute(bu - bl,[1 3 2]);
            % Map bounds of the factors to bounds of the constraint 
            % zonotope.
            c_ = c + reshape(pagemtimes(G,bc),[n bSz]);
            r = dr + reshape(pagemtimes(abs(G),br),[n bSz]);
            l = c_ - r;
            u = c_ + r;

            % Bound the solution of the linear program.
            rl = -aux_boundConZonotopeSupportFunction(A,b,-G,options);
            % Bound the solution of the linear program.
            ru = aux_boundConZonotopeSupportFunction(A,b,G,options);

            % Compute the bounds of the constrained zonotope from the 
            % bounded solutions.
            l = max(l,c + rl - dr);
            u = min(u,c + ru + dr);
        case 'exact'
            % Slow implementation that computes the exact bounds by solving 
            % 2*n linear programs.
    
            % Initialize the bounds.
            l = NaN([n bSz]);
            u = NaN([n bSz]);
            bl = NaN([n bSz]);
            bu = NaN([n bSz]);
    
            % Iterate over the batch entries.
            for i=1:bSz
                % Obtain parameters of the i-th batch entry.
                Ai = double(gather(A(:,:,i)));
                bi = double(gather(b(:,i)));
                ci = double(gather(c(:,i)));
                Gi = double(gather(G(:,:,i)));
                if all(~isnan(Ai),'all') && all(~isnan(bi),'all')
                    % Construct linear program.
                    prob = struct('Aineq',Ai,'bineq',bi, ...
                        'lb',-ones([q 1]),'ub',ones([q 1]));
                    % Loop over the dimensions.
                    for j=1:n
                        % Find the lower bound for the j-th dimension.
                        prob.f = Gi(j,:);
                        % Solve the linear program.
                        [~,lij,efl] = CORAlinprog(prob);
                        % Find the upper bound for the j-th dimension.
                        prob.f = -Gi(j,:);
                        % Solve the linear program
                        [~,uij,efu] = CORAlinprog(prob);

                        if efl > 0 && efu > 0
                            % Solutions found; assign values.
                            l(j,i) = ci(j) + lij;
                            u(j,i) = ci(j) - uij;
                        else
                            % No solution; the constrained zonotope is empty.
                            break;
                        end
                    end

                    % Loop over the generators.
                    for j=1:q
                        % Find the lower bound for the j-th dimension of 
                        % the constrained hypercube.
                        prob.f = double((1:q) == j);
                        % Solve the linear program.
                        [~,blij,efbl] = CORAlinprog(prob);
                        % Find the upper bound for the j-th dimension of 
                        % the constrained hypercube.
                        prob.f = -double((1:q) == j);
                        % Solve the linear program
                        [~,buji,efbu] = CORAlinprog(prob);

                        if efbl > 0 && efbu > 0
                            % Solutions found; assign values.
                            bl(j,i) = blij;
                            bu(j,i) = -buji;
                        else
                            % No solution; the constrained zonotope is empty.
                            break;
                        end
                    end
                end
            end
            % Add the approximation error.
            l = l - dr;
            u = u + dr;
        otherwise 
            % Invalid option.
            throw(CORAerror('CORA:wrongFieldValue', ...
                'options.nn.conzonotope_bounding_method', ...
                    {'fourier-motzkin','dual-iter','exact'}));
    end
end

function [bl,bu] = aux_fourierMotzkinApproximation(A,b,options)
    % Compute the bounds [bl,bu] of a bounded polytope P:
    % Given P=(A,b) \cap [-1,1], compute its bounds, i.e., 
    % [bl,bu]\supseteq {x\in\R^q\mid A\,x\leq b} \cap [-1,1].

    % Specify a numerical tolerance to avoid numerical instability.
    tol = 1e-8;

    % Efficient approximation by isolating the i-th variable. -------------
    % We compute a box-approximation of the valid factor for the 
    % constraint zonotope, 
    % i.e., [\underline{\beta},\overline{\beta}] 
    %   \supseteq \{\beta \in [-1,1]^q \mid A\,\beta\leq b\}.
    % We view each constraint separately and use the tightest 
    % bounds. For each constraint A_{(i,\cdot)}\,\beta\leq b_{(i)}, 
    % we isolate each factor \beta_{(j)} and extract bounds:
    % A_{(i,\cdot)}\,\beta\leq b_{(i)} 
    %   \implies A_{(i,j)}\,\beta_{(j)} \leq 
    %       b_{(i)} - \sum_{k=1,...,q, k\neq j} A_{(i,k)}\,\beta_{(k)}
    % Based on the sign of A_{(i,j)} we can either tighten the 
    % lower or upper bound of \beta_{(j)}.
    
    % Specify maximum number of iterations.
    maxIter = options.nn.polytope_bound_approx_max_iter;
    
    % Initialize bounds of the factors.
    bl = -ones(size(A,[2 3]),'like',A);
    bu = ones(size(A,[2 3]),'like',A);
    
    % Permute the dimension of the constraints for easier handling.
    A_ = permute(A,[2 1 3]);
    b_ = permute(b,[3 1 2]);
    % Reshape factor bounds for easier multiplication.
    bl_ = permute(bl,[1 3 2]);
    bu_ = permute(bu,[1 3 2]);
        
    % Initialize iteration counter.
    iter = 1;
    tighterBnds = 1;
    while tighterBnds && iter <= maxIter
        % Scale the matrix entries with the current bounds.
        ABnd = A_.*((A_ > 0).*bl_ + (A_ < 0).*bu_);
        % Isolate the i-th variable of the j-th constraint. Sum all
        % variables execpt the i-th. Compute right-hand side of the 
        % inequalities.
        rh = min(max((b_ - (sum(ABnd,1) - ABnd))./A_,bl_),bu_);
        % Update the lower bounds.
        bnd = repmat(bl_,1,size(A,1),1);
        bnd(A_ < 0) = rh(A_ < 0);
        bl_ = max(bnd,[],2);
        % Update the upper bounds.
        bnd = repmat(bu_,1,size(A,1),1);
        bnd(A_ > 0) = rh(A_ > 0);
        bu_ = min(bnd,[],2);
        % Check if the bounds could be tightened.
        tighterBnds = any( ...
            (bl + tol < bl_(:,:) | bu_(:,:) < bu - tol) ... tighter bounds
                & bl_(:,:) <= bu_(:,:), ... not empty
            'all');
        bl = reshape(bl_, size(bl));
        bu = reshape(bu_, size(bu));
        % Increment iteration counter.
        iter = iter + 1;
    end
end

function h = aux_boundConZonotopeSupportFunction(A,b,G,options)
    % We compute bounds of a constrained zonotope by computing a bound on 
    % the support function in both directions for each dimension. 
    % For a constrained zonotope (c,G,A,b) we the bounds 
    % [l_{(i)},u_{(i)}] for the i-th dimension are
    %
    % [Lower Bound]     l_{(i)} = c_{(i)} + min G_{(i,.)}\,\beta 
    %                                       s.t. A\,\beta <= b
    %                                              -\beta <= 1
    %                                               \beta <= 1
    %
    % [Upper Bound]     u_{(i)} = c_{(i)} - min -G_{(i,.)}\,\beta 
    %                                       s.t. A\,\beta <= b
    %                                              -\beta <= 1
    %                                               \beta <= 1
    %
    % We bound support function of a constrained zonotope by optimizing the
    % dual variables. By the weak duality every dual solution is 
    % automatically a bound on the primal solution. The dual problems with 
    % variables \lambda are:
    % (Additional variables \mu and \nu are required for the bounds of
    % \beta, i.e., -1 <= \beta <= 1.)
    %
    % [Lower Bound]     max -b^T\,\lambda - 1^T\,(\mu + \nu)
    %                   s.t. A^T\,\lambda - \mu + \nu = -G_{(i,.)}^T
    %                               \lambda >= 0, \mu >= 0, \nu >= 0
    % 
    % [Upper Bound]     max -b^T\,\lambda - 1^T\,(\mu + \nu)
    %                   s.t. A^T\,\lambda - \mu + \nu = G_{(i,.)}^T
    %                               \lambda >= 0, \mu >= 0, \nu >= 0
    %
    % We can reduce the variables \mu and \nu: For a fixed \lambda >= 0; 
    % from the constraints, we have
    %
    % [Lower Bound] -\mu + \nu = -A^T\,\lambda - G_{(i,.)}^T, and 
    %
    % [Upper Bound] -\mu + \nu = -A^T\,\lambda + G_{(i,.)}^T. 
    %
    % Therefore, \mu and \nu are optimized by (for both bounds)
    %
    % [Lower Bound]    max -1^T\,(\mu + \nu) = -1^T\,\abs(-A^T\,\lambda - G_{(i,.)}^T)\,1
    %                  s.t. \mu,\nu >= 0                
    % 
    % [Upper Bound]    max -1^T\,(\mu + \nu) = -1^T\,\abs(-A^T\,\lambda + G_{(i,.)}^T)\,1
    %                  s.t. \mu,\nu >= 0
    %
    % Thus, we can use the weak duality to obtain bounds on l_{(i)} and 
    % u_{(i)} with any \lambda >= 0 by
    %
    % [Lower Bound]    l_{(i)} >= c_{(i}) - min_{\lambda >= 0} \lambda^T\,b + \abs{-G_{(i,.)}^T - A^T\,\lambda}\,1 
    %
    % [Upper Bound]    u_{(i)} <= c_{(i}) + min_{\lambda >= 0} \lambda^T\,b + \abs{G_{(i,.)}^T - A^T\,\lambda}\,1  
    %
    % We use gradient-based optimization to optimize \lambda. The problem
    % is convex, but non-smooth, i.e., the optimization objective is
    %   \phi(\lambda) = \lambda^T\,b + \abs{+/-G_{(i,.)}^T - A^T\,\lambda}\,1,
    % for \lambda >= 0. 
    % We approximate the objective with a smooth approximation of \abs{.}, 
    % i.e., 
    %   \abs{x} \approx \sqrt{x^2 + epsilon^2} - epsilon
    % The surrogate objective is
    %   \tilde{\phi}(\lambda) = \lambda^T\,b + (\sqrt{(+/-G_{(i,.)}^T - A^T\,\lambda)^2 + epsilon^2} - epsilon)\,1,
    % The gradient of \tilde{\phi} is
    %   \nabla_{\lambda} \tilde{\phi}(\lambda) = b - A\,/\sqrt{+/-G_{(i,.)}^T - A^T\,\lambda + epsilon^2}\,1.

    % Obtain the number of dimensions, generators, and the batch size.
    [n,q,bSz] = size(G);
    % Obtain the number of constraints.
    [p,~,~] = size(A);
    % Reshape for easier batch-wise computations.
    b = reshape(b,[p 1 bSz]);

    % Compute scaling factors for the constraints (rows of A).
    % This prevents certain constraints from dominating the gradient 
    % purely due to magnitude.
    scalingA = 1./(max(1,sqrt(sum(A.^2,2)))); 
    A = scalingA.*A;
    b = scalingA.*b;

    % Specify an optimization step size.
    sSz = options.nn.conzonotope_bound_step_size;
    % Specify step size decay.
    beta = 0.5;
    % Specify the number of step sizes to try during backtracking.
    nsSzs = 2^5;
    % Compute the possible step sizes. Permute s.t. we can do batch wise
    % computations for all step sizes simultaneously.
    sSzs = permute((beta.^(0:(nsSzs-1)))*sSz,[1 3 4 2]);

    % Obtain the number of iterations.
    maxIter = options.nn.conzonotope_bound_max_iter;

    % Specify a smoothing parameter.
    epsilon = 1;
    % Construct a function handle for the surrogate absolute value.
    % softplus = @(x) log(1 + exp(x));
    % dsoftplus = @(x) 1./(1 + exp(-x));
    sAbs = @(x,epsilon) sum(sqrt(x.^2 + epsilon^2) - epsilon,2);
    dsAbs = @(x,epsilon) x./sqrt(x.^2 + epsilon^2);
    % Define strong convexity parameter (1e-6 to 1e-4 is usually safe)
    mu = 1e-5;
    % Construct a function handle for the optimization objective.
    f = @(x) reshape(pagemtimes(x,'transpose',b,'none') ...
        + sum(abs(G - pagemtimes(x,'transpose',A,'none')),2),n,bSz,[]) ...
        + (mu/2)*reshape(sum(x.^2,1),n,bSz,[]);
    df = @(x,epsilon) b - pagemtimes(A,'none', ...
        dsAbs(G - pagemtimes(x,'transpose',A,'none'),epsilon),'transpose') ...
        + mu*x;

    % Initialize the dual variables.
    l = zeros([p n bSz],'like',G); % zeros([p n bSz],'like',G);
    % Initialize iteration counter.
    iter = 1;

    % Initialize the L-BFGS parameters and memory.
    mSz = 5; % memory size.
    c1 = 1e-4; % for backtracking line search (Armijo condition).
    c2 = 0.9; % for backtracking line search (Curvature condition).
    % Initialize the memory; we store in reverse order, i.e., 
    % k-1, k-2, ..., k-mSz.
    s = nan([p n bSz 0],'like',G); % store variables differences
    y = nan([p n bSz 0],'like',G); % store gradient differences
    % Store intermediate computed factors.
    as = nan([1 n bSz mSz],'like',G);
    bs = nan([1 n bSz mSz],'like',G);

    % Compute the gradient at the initial position.
    dl = df(l,epsilon);

    % Do an iterative optimization of the dual variables.
    while iter <= maxIter
        % Store old variables and gradients.
        l_ = l;
        dl_ = dl;

        % Estimate the search direction using two-recursions.
        qk = dl;
        for i=1:min(mSz,size(s,4))
            ri = 1./sum(y(:,:,:,i).*s(:,:,:,i),1);
            as(:,:,:,i) = ri.*sum(s(:,:,:,i).*qk,1);
            qk = qk - as(:,:,:,i).*y(:,:,:,i);
        end
        if isempty(s)
            % The memory is empty.
            gammak = 1;
        else
            gammak = sum(s(:,:,:,1).*y(:,:,:,1),1)./...
                sum(y(:,:,:,1).*y(:,:,:,1),1);
        end
        zk = gammak.*qk; % Hk = gammak.*eye(p,'like',G);
        for i=min(mSz,size(s,4)):-1:1
            ri = 1./sum(y(:,:,:,i).*s(:,:,:,i),1);
            bs(:,:,:,i) = ri.*sum(y(:,:,:,i).*zk,1);
            zk = zk + s(:,:,:,i).*(as(:,:,:,i) - bs(:,:,:,i));
        end
        % Invert final search direction for minimization.
        zk = -zk;

        % Zero search direction if l == 0 and dl > 0.
        zk(l == 0 & zk < 0) = 0;

        % Do a backtracking line search to find a suitable step size. 
        % Therefore, we do a simultaneous step with all step sizes and 
        % pick the best result.

        % Do a proximal gradient step with the current step size.
        ls = l + sSzs.*zk;
        % Project onto the feasible set, i.e., lambda >= 0.
        ls = max(ls,0);

        % Compute objective at the old and new point.
        fl = f(l);
        fls = f(ls);
        dls = df(ls,epsilon);
        % Compute the scalar product between gradient and search direction.
        sdflzk = reshape(sum(dl.*zk,1),[n bSz]);
        sdflszk = reshape(sum(dls.*zk,1),[n bSz nsSzs]);
        % Compute the Armijo condition.
        armijoCond = (fls <= fl + c1.*sSzs(:,:,:).*sdflzk);
        curvCond = (abs(sdflszk) <= c2.*abs(sdflzk));

        % Use 'max' function to find the first, i.e., largest, step size
        % that sastisfies Armijo–Goldstein the condition.
        [~,ids] = max(armijoCond & curvCond,[],3);
        % Compute linear indices to extract the correct values.
        idx = sub2ind([n bSz nsSzs], ...
            repelem((1:n)',1,bSz),repmat(1:bSz,n,1),ids);
        % Update the variables with the largest possible step size.
        l = reshape(ls(:,idx),[p n bSz]);

        % Compute the gradient at the new position; projected gradient.
        dl = df(l,epsilon); % (l - l_)./sSzs(ids); % 

        if all(curvCond,'all')
            % Only update the memory if the curvature condition holds.

            % Compute new memory entries.
            s_ = l - l_;
            y_ = dl - dl_;
            % Prepend the new entries to the memory.
            s = cat(4,s_,s(:,:,:,1:min(mSz-1,end)));
            y = cat(4,y_,y(:,:,:,1:min(mSz-1,end)));
        end

        % Increment iteration counter.
        iter = iter + 1;
        % Decay epsilon.
        epsilon = max(0.5*epsilon,1e-6);
    end
    % Compute the final bound on the support function.
    h = f(l);
end

% ------------------------------ END OF CODE ------------------------------
