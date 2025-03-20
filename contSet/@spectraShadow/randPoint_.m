function p = randPoint_(SpS,N,type,varargin)
% randPoint_ - generates a random point within a spectrahedral shadow
%
% Syntax:
%    p = randPoint_(SpS,N,type)
%
% Inputs:
%    SpS - spectraShadow object
%    N - number of random points
%    type - type of the random point ('standard' or 'uniform')
%
% Outputs:
%    p - array of random points in SpS
%
% Example: 
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%    points = randPoint(SpS,100);
%
%    figure; hold on;
%    plot(SpS,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, polytope/randPoint_

% Authors:       Adrian Kulmburg, Vladimir Popa
% Written:       12-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check if empty
    if representsa(SpS,'emptySet',1e-5)
        p = double.empty(dim(SpS),0); return
    end
    
    if ~isBounded(SpS)
        % Currently, the methods implemented only work for bounded sets
        throw(CORAerror('CORA:specialError',...
            ['The function randPoint is currently only implemented for '...
            'bounded spectrahedra/spectrahedral shadows.'])); 
    end

    % Special case: SpS is just a point
    [isPoint, p] = representsa(SpS, 'point', 1e-5);
    if isPoint
        p = repmat(p, N);
        return
    end
    
    % Choose the appropriate method
    if strcmp(type,'standard') || strcmp(type,'uniform') ...
            || strcmp(type,'uniform:hitAndRun')
        % The default is hit-and-run, since it gave the best results so far
        p = aux_randPointHitAndRun(SpS, N);
    elseif strcmp(type,'extreme')
        p = aux_randPointExtreme(SpS, N);
    elseif strcmp(type, 'uniform:ballWalk')
        p = aux_randPointBallWalk(SpS, N);
    elseif strcmp(type, 'uniform:billiardWalk')
        p = aux_BilliardWalk(SpS, N);
    else
        throw(CORAerror('CORA:noSpecificAlg',type,SpS));
    end
end


% Auxiliary functions -----------------------------------------------------

function p = aux_BilliardWalk(SpS, N)
    % Allocate Memory
    p = zeros(dim(SpS),N);
    
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
    end

    G = SpS.G;
    c = SpS.c;
    [A0,Ai] = priv_getCoeffMatrices(SpS);
    
    % Need to check whether the spectrahedral shadow is degenerate
    [res, X] = isFullDim(SpS);
    r = dim(SpS);
    if ~res
        X = [X zeros([r r-size(X,2)])];
    end

    [~, p0] = supportFunc_(SpS - c, rand([dim(SpS) 1]), 'upper'); % out starting point
    R0 = 10 * dim(SpS); % reflection number
    n = size(G, 1);
    m = size(G, 2);

    for i=1:N
        
        % segment start point at p0
        q0 = p0;
        % length l
        I = interval(SpS);
        tau = max(supremum(I)-infimum(I));
        l = tau * -log(rand(1));
        % reflection number
        R = R0;

        % Select a random direction
        direction = randn([dim(SpS) 1]);
        direction = direction./norm(direction);
        
        % in case the spectral shadow is degenerate, the direction can only
        % be in the subspace spanned by the vectors of X (which are
        % normalized and orthogonal, so all is well, we still have a
        % uniformly sampled direction)
        direction = X * direction;

        while true
            
            % Set up optimization problem to compute the maximal value of
            % lambda (i.e., how far we can go in direction, starting from
            % q0)
            lambda = sdpvar(1, 1);
            x = sdpvar(m, 1);

            constraints = [q0 + lambda*direction == G*x];

            lmi = A0;
            for j=1:m
                lmi = lmi + Ai{j}*x(j);
            end

            constraints = [constraints, lmi>=0];

            % Solve optimization problem
            diagnostics = optimize(constraints, -lambda, options);

            lambda = value(lambda);

            y = q0 + lambda * direction;

            % compute the normal vector maxX only if 
            if l > lambda
                % find a boundary point y in SpS
                % argmax y^Tp s.t. y is in SpS*
                % min_x,Z y^Tx s.t. Tr(A0 * Z) <= 1, Tr(Ai * Z) = ei^T * G^T * x, Z >= 0
                x = sdpvar(n, 1);
                Z = sdpvar(size(A0, 1), size(A0, 1));
        
                constraints = [trace(A0' * Z) <= 1];
    
                % canonnical base vectors
                e=@(k,n) [zeros(k-1,1);1;zeros(n-k,1)];
         
                for j=1:m
                    constraints = [constraints, trace(Ai{j} * Z) == e(j,m)' * G' *x];
                end
    
                constraints = [constraints, Z >= 0];
    
                objective = x' * y;
                
                diagnostics = optimize(constraints, objective,options);
                
                % sol{1} is optimal x
                % sol{2} is optimal Z
                maxX = value(x);
                
                s = maxX / norm(maxX);
                % change start point, direction and remaining length
                q0 = y;
                direction = direction - 2 * dot(direction,s) * s;
                l = l - lambda;
            else
                p0 = q0 + direction * l;
                p(:,i) = p0;
                break;
            end

            % sample new direction from the last returned point if the
            % maximum reflection number
            if R <= 0
                R = R0;
                q0 = p0;
                l = tau *-log(rand(1));
                direction = randn([dim(SpS) 1]);
                direction = direction./norm(direction);

                direction = X * direction;
                continue;
            else
                R = R - 1;
            end
        end
    end
    p = p + c;
end

function p = aux_randPointBallWalk(SpS, N)
    % Allocate Memory
    p = zeros(dim(SpS), N);

    % Need to check whether the spectrahedral shadow is degenerate
    [res, X] = isFullDim(SpS);
    r = dim(SpS);
    if ~res
        X = [X zeros([r r-size(X,2)])];
    end

    I = interval(SpS);

    l = (I.supremum()-I.infimum())./2;
    c = center(I);
    SpS_0 = SpS - c;
    SpS0_trans = diag(l.^(-1))*SpS_0;

    % Radius of the Ball
    r = 1/sqrt(dim(SpS));

    [~, p0] = supportFunc_(SpS0_trans,2.*rand([dim(SpS) 1])-1,'upper');

    y = p0;
    for i=1:N

        doWhile = false; % Need to pass the next loop at least once
        while ~doWhile || ~SpS0_trans.contains(y,'exact',1e-8)
            doWhile = true;

            % Generate a vector of dim random numbers
            direction = randn([dim(SpS0_trans) 1]);
            direction = direction./norm(direction);

            % in case the spectrahedral shadow is degenerate, the direction can only be in
            % the subspace spanned by the vectors of X (which are normalized
            % and orthogonal, so all is well, we still have a uniformly sampled
            % direction)
            direction = X * direction;

            % Generate a random number
            d = rand(1);

            % Benerate x in the Ball of radius r
            x = d^(1/dim(SpS)) * direction;

            y = p0 + r * x;

                
        end

        p(:,i) = y;
        p0 = y;
    end

    p = diag(l)*p;
    p = p + c;
   
end

function p = aux_randPointHitAndRun(SpS, N)
    % Allocate Memory
    p = zeros(dim(SpS),N);
    
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
    end
    
    G = SpS.G;
    c = SpS.c;
    [A0,Ai] = priv_getCoeffMatrices(SpS);
    
    % Need to check whether the polytope is degenerate
    [res, X] = isFullDim(SpS);
    n = dim(SpS);
    if ~res
        X = [X zeros([n n-size(X,2)])];
    end
    
    % We start at some random point in the spectrahedral shadow; for this, 
    % we can use the support function in a random direction:
    [~,currentPoint] = supportFunc_(SpS,rand([dim(SpS) 1]),'upper');
    for i=1:N
        % Select a random direction
        direction = randn([dim(SpS) 1]);
        direction = direction./norm(direction);
        
        % in case the spectral shadow is degenerate, the direction can only
        % be in the subspace spanned by the vectors of X (which are
        % normalized and orthogonal, so all is well, we still have a
        % uniformly sampled direction)
        direction = X * direction;
        
        % Compute the range of the distance we can go in the selected
        % direction, without leaving the spectrahedral shadow. This means
        % that we are searching for the allowed lambda, such that
        % currentPoint+lambda*direction \in S
        lambda = sdpvar(1,1);
        beta = sdpvar(size(G,2),1);

        % residual to avoid numerical errors; not that this, in turn, can
        % lead to points that are not quite exactly within the
        % spectrahedral shadow, but they are up to a certain tolerance
        s = sdpvar(size(A0,1), 1);
        residual = diag(s);
        tol = 1e-7;
        constraints_residual = [];
        for j=1:size(A0,1)
            constraints_residual = [constraints_residual s(j)<=tol -s(j)<=tol];
        end
        
        % Set up the rest of the optimization problem
        constraints = G * beta + c == currentPoint + lambda * direction;
        LMI = A0;
        for j=1:size(G,2)
            LMI = LMI + beta(j) * Ai{j};
        end
        constraints = [constraints LMI + residual>=0];
        constraints = [constraints constraints_residual];
        
        cost = lambda;

        diagnostics = optimize(constraints,cost,options);
        exitflag = diagnostics.problem;
        if exitflag == 1
            % If the problem is infeasible, this means that the solver,
            % despite our best efforts, wasn't smart enough to see that the
            % solution lambda = 0 was available, so we need to set it
            % ourselves
            lambda_lower = 0;
        else
            lambda_lower = value(lambda);
        end

        cost = -lambda;

        diagnostics = optimize(constraints,cost,options);
        exitflag = diagnostics.problem;
        
        if exitflag == 1
            % If the problem is infeasible, this means that the solver,
            % despite our best efforts, wasn't smart enough to see that the
            % solution lambda = 0 was available, so we need to set it
            % ourselves
            lambda_upper = 0;
        else
            lambda_upper = value(lambda);
        end


        % Select random distance:
        lambda = rand(1)*(lambda_upper-lambda_lower) + lambda_lower;
        %lambda = unifrnd(lambda_lower, lambda_upper);
        
        % Compute next point:
        currentPoint = currentPoint + lambda * direction;
        
        % Add it to the list:
        p(:,i) = currentPoint;
    end
end

function p = aux_randPointExtreme(SpS, N)
    % Start with sampling some uniform points
    p = aux_randPointHitAndRun(SpS,N);
    
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
    end
    
    G = SpS.G;
    c = SpS.c;
    [A0,Ai] = priv_getCoeffMatrices(SpS);
    
    % Need to check whether the polytope is degenerate
    [res, X] = isFullDim(SpS);
    n = dim(SpS);
    if ~res
        X = [X zeros([n n-size(X,2)])];
    end

    % For each point, choose a random direction, and go as far in that
    % direction as possible, starting from p
    for i=1:N
        
        % Select a random direction
        direction = randn([dim(SpS) 1]);
        direction = direction./norm(direction);
        
        % in case the spectral shadow is degenerate, the direction can only
        % be in the subspace spanned by the vectors of X (which are
        % normalized and orthogonal, so all is well, we still have a
        % uniformly sampled direction)
        direction = X * direction;
        
        % Compute the range of the distance we can go in the selected
        % direction, without leaving the spectrahedral shadow. This means
        % that we are searching for the allowed lambda, such that
        % currentPoint+lambda*direction \in S
        lambda = sdpvar(1,1);
        beta = sdpvar(size(G,2),1);
        
        constraints = G * beta + c == p(:,i) + lambda * direction;
        LMI = A0;
        for j=1:size(G,2)
            LMI = LMI + beta(j) * Ai{j};
        end
        constraints = [constraints LMI>=0];
        
        cost = -lambda;
        diagnostics = optimize(constraints,cost,options);
        
        lambda = value(lambda);
     
        % Compute next point:
        p(:,i) = p(:,i) + lambda * direction;
        
    end
end

% ------------------------------ END OF CODE ------------------------------
