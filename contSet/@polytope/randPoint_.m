function pt = randPoint_(P,N,type,varargin)
% randPoint_ - generates a random point within a polytope
%
% Syntax:
%    pt = randPoint_(P)
%    pt = randPoint_(P,N)
%    pt = randPoint_(P,N,type)
%    pt = randPoint_(P,'all','extreme')
%
% Inputs:
%    P - polytope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    pt - random point in R^n
%
% Example: 
%    P = polytope.generateRandom('Dimension',2);
%    points = randPoint(P,100);
%
%    figure; hold on;
%    plot(P,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, zonotope/randPoint_

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Adrian Kulmburg
% Written:       30-October-2020
% Last update:   25-June-2021 (MP, add type gaussian)
%                25-May-2023 (AK, added uniform sampling)
%                18-January-2024 (MW, fix unbounded case in extreme)
% Last revision: 04-April-2023 (MW, restructure)

% ------------------------------ BEGIN CODE -------------------------------

    % fullspace
    if representsa_(P,'fullspace',0)
        % sample according to normal distribution
        pt = randn(dim(P),N); return
    end

    % empty
    if representsa(P,'emptySet')
        pt = zeros(dim(P),0); return
    end
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        pt = vertices(P); return
    end
    
    % init random points
    pt = zeros(dim(P),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            pt(:,i) = aux_randPointStandard(P);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            pt(:,i) = aux_randPointExtreme(P);
        end
        
    % Default algorithm for the uniform distribution is the billiard walk,
    % since it gives more accurate results, despite a slightly longer
    % computation time
    elseif strcmp(type,'uniform') || strcmp(type,'uniform:billiardWalk')
        pt = aux_randPointBilliard(P, N);
    elseif strcmp(type,'uniform:hitAndRun')
        pt = aux_randPointHitAndRun(P, N);
    else
        throw(CORAerror('CORA:noSpecificAlg',type,P));
    end
end


% Auxiliary functions -----------------------------------------------------

function pt = aux_randPointStandard(P)
% generate random point within the polytope

    % draw n+1 random extreme points
    n = dim(P);
    points = zeros(n,n+1);
    
    for i = 1:size(points,2)
        points(:,i) = aux_randPointExtreme(P); 
    end

    % interpolate between the points
    delta = rand(size(points,2),1);
    delta = delta./sum(delta);
    
    pt = points*delta;
    
end

function pt = aux_randPointExtreme(P)
% generate random point on the polytope boundary

    if ~isBounded(P)
        % go along halfspace vectors to find finite-valued extreme points
        A = [P.A; P.Ae; -P.Ae]; b = [P.b; P.be; -P.be];
        picked_halfspace = randi([1,length(b)]);
        [~,pt] = supportFunc_(P,A(picked_halfspace,:),'upper');
        
    else
        % select random direction, normalize
        n = dim(P);
        d = rand(n,1) - 0.5*ones(n,1);
        d = d./norm(d);

        % compute farthest point in this direction still inside the polytope
        [~,pt] = supportFunc_(P,d,'upper');
    end
    
end

function p = aux_randPointBilliard(P, N)
    
    % Heuristic choices for tau and R
    tau = norm(interval(P),2,'ub');
    R = 10 * dim(P);


    %polytope initialization

    p = zeros(dim(P),N);
    A = P.A;
    b = P.b;
    
    % Need to check whether the polytope is degenerate
    n = dim(P);
    [res, X] = isFullDim(P);
    if ~res
        X = [X zeros([n n-size(X,2)])];
    end

    %starting random point
    startingPoint = randPoint(P);

    currentPoint = startingPoint;
    
    % This will be the amount of points we have generated so far
    i = 0;
    
    while i<N % while we have generated fewer than N points
        
        % generate a sampled length, using the formula given in the paper
        ell = -tau * log(rand(1));
        
        %generate a random direction 'd'
        direction = randn([dim(P) 1]);
        direction = direction./norm(direction);
        
        % in case the polytope is degenerate, the direction can only be in
        % the subspace spanned by the vectors of X (which are normalized
        % and orthogonal, so all is well, we still have a uniformly sampled
        % direction)
        direction = X * direction;
        
        % set the number of remaining reflections to R
        remainingReflections = R;

        % Once we have selected a total length (that's this ell above),
        % and selected a direction it should go to,
        % we can try to construct the billiard path. This billiard path
        % will have a final destination, which we call like this:
        billiard_path = currentPoint;
        % At the start, the constructed destination is just the
        % currentPoint (i.e., the previous point we ended up with).
        % But right now, we have a path length of zero:
        path_length = 0;
        % Our goal is to compute the billiard-path with length = ell,
        % starting off in direction 'direction'
        while path_length < ell && remainingReflections > 0
            % This is now going to be the bulk of the algorithm. Here, we
            % do what you are familiar with, so we compute the lambdas, etc...

            %generate a lambdas 'l'
            lambdas = (b-A*billiard_path)./(A*direction);
            lambdas = lambdas(~isnan(lambdas));

            %p point
            lambda_upper = min(lambdas(A*direction >= 0));

            % (We only need lambda_upper)
            
            % In the case where the remaining path length is small
            % enough...
            remaining_length = ell - path_length;
            remainingReflections = remainingReflections - 1;
            if remaining_length <= lambda_upper
                % ... then we do not need reflections, and we can just
                % find the next destination with the remaining path
                % length...
                billiard_path = billiard_path + remaining_length .* direction;
                % ... say that we found a path of proper length ...
                path_length = ell;
                % ... update the current iteration ...
                currentPoint = billiard_path;
                % ... add it to our output ...
                p(:,i+1) = billiard_path;
                % ... and say that we have found a new point.
                i = i+1;
                % We can then just break of the while-loop, because we
                % were successful.
                break
            end
            
            % If the remaining length exceeds lambda_upper, it means that
            % we would need to go beyond the boundary, which is not
            % allowed. In that case, we need to use reflections.

            % To make a reflection, we need to know what facet of the
            % polytope will do the reflection. To do that, we need to
            % check which half-space is reached 'first', and so is
            % responsible for lamda_upper. There is probably a smarter
            % way to do this, but for now consider this: Since we know
            % lambda_upper, we know that the point on the boundary is
            % given as...
            p_boundary = billiard_path + lambda_upper .* direction;
            
            % Now, we need to check on which half-space p_boundary is.
            % The best way to do this is to plug it into the H-rep, and
            % check which inequality becomes an equality:
            index = find(round(A*p_boundary-b, 8) == 0);
            % The rounding is just to make things work well enough; if
            % you don't put it, some floating point errors might appear.
            
            % Ideally, p_boundary is exactly on one half-space. If that
            % is not the case, it means that it is on an edge of the
            % polytope, which has probability zero. Still, it might
            % happen, and if it does we need to cancel the current 
            % iteration and start again. The easiest way to do this is to
            % set path_length to infinity, or you could also break:
            if length(index) > 1
                path_length = inf;
                break
            end
            
            % So, assuming all went well and we only have one half-space,
            % using index we know what half-space we need:
            halfspace_orthogonal = A(index,:)'./norm(A(index,:));
            % The rows of A give the vector that is orthogonal to the
            % half-space
            
            % Now, we can do the reflection itself. To do that, we first
            % update the current billiard_path, to the point on the
            % boundary:
            billiard_path = p_boundary;
            % By doing that, we went through a line of length
            % lambda_upper, so the current path_length also has to be
            % updated:
            path_length = path_length + lambda_upper;
            
            % Finally, we need to update the direction. This is given
            % using the formula of the paper:
            direction = direction - 2 .* (direction' * halfspace_orthogonal) .* halfspace_orthogonal;
            % Just for good measure, we normalize it, to make sure
            % nothing bad can happen:
            direction = direction./norm(direction);
            
            % And that's it, we have made a reflection :) Since we did
            % not add any points to the output, the while-loop will be used
            % at least once again to continue the reflection, until it
            % either hits an edge or finishes the path.
        end
        
    end


end

function p = aux_randPointHitAndRun(P, N)
    % Allocate Memory
    p = zeros(dim(P),N);
    
    A = P.A;
    b = P.b;
    % This is just a convenient renaming. The polytope is now the set of
    % points x, s.t. Ax <= b.
    
    % Need to check whether the polytope is degenerate
    [res, X] = isFullDim(P);
    n = dim(P);
    if ~res
        X = [X zeros([n n-size(X,2)])];
    end
    
    % We start at some random point in the polytope; for this, we can use
    % the other, easier method:
    currentPoint = randPoint(P);
    for i=1:N
        % Select a random direction
        direction = randn([dim(P) 1]);
        direction = direction./norm(direction);
        
        % in case the polytope is degenerate, the direction can only be in
        % the subspace spanned by the vectors of X (which are normalized
        % and orthogonal, so all is well, we still have a uniformly sampled
        % direction)
        direction = X * direction;
        
        % Compute the range of the distance we can go in the selected
        % direction, without leaving the polytope. This means that we are
        % searching for the allowed lambda, such that
        % A(currentPoint+lambda*direction) <= b
        % Doing this componentwise, we get for all possible edge-cases for
        % lambda:
        lambdas = (b-A*currentPoint)./(A*direction);
        
        % Remove all 'useless' values
        lambdas = lambdas(~isnan(lambdas));
        
        % Now, the lambdas that are 'upper bounds' correspond to the cases,
        % where (A*direction)(i) is positive (it doesn't change the sign of
        % the inequality):
        lambda_upper = min(lambdas(A*direction >= 0));
        
        % Similarly, the lambdas that are 'lower bounds' correspond to the
        % cases, where (A*direction)(i) is negative (it changes the sign):
        lambda_lower = max(lambdas(A*direction < 0));
        
        % Checking that the polytope is well defined:
        if abs(lambda_upper) == inf || abs(lambda_lower) == inf
            throw(CORAerror('CORA:specialError','The polytope is either empty or unbounded, therefore one cannot sample it uniformly.'));
        end
        % Checking that lambda_upper >= lambda_lower
        if lambda_upper < lambda_lower
            if abs(lambda_upper-lambda_lower) < 1e-6
                % It could be that currentPoint is on an edge/a vertex, and
                % then it could be that lambda_upper and lambda_lower are
                % wrongly calculated due to floating point errors. If that
                % is the case, just choose lambda_upper=lambda_lower=0
                lambda_upper = 0;
                lambda_lower = 0;
            else
                throw(CORAerror('CORA:specialError','Something unexpected happened during the execution of the hit-and-run algorithm.'));
            end
        end
        
        % Select random distance:
        lambda = unifrnd(lambda_lower, lambda_upper);
        
        % Compute next point:
        currentPoint = currentPoint + lambda * direction;
        
        % Add it to the list:
        p(:,i) = currentPoint;
    end
end


% ------------------------------ END OF CODE ------------------------------
