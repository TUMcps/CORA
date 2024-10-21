function res = contains_(SpS,obj,type,tol,maxEval,varargin)
% contains_ - determines if a spectrahedral shadow contains a set or a
%    point
%
% Syntax:
%    res = contains_(SpS,obj)
%
% Inputs:
%    SpS - spectraShadow object
%    obj - contSet object or vector
%    type - not implemented yet
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example:
%    A0 = eye(4);
%    A1 = [-1 0 0 0;0 1 0 0;0 0 0 0;0 0 0 0];
%    A2 = [0 0 0 0;0 0 0 0;0 0 -1 0;0 0 0 1];
%    SpS = spectraShadow([A0 A1 A2]);
%
%    SpS.contains([0;0]) % true
%    SpS.contains([2;2]) % false
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains_, conZonotope/contains_

% Authors:       Maximilian Perschl, Adrian Kulmburg
% Written:       08-June-2023 
% Last update:   ---             
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


    if isnumeric(obj)
        res = aux_containsPoint(SpS, obj, tol);
    else
        try
            obj = spectraShadow(obj);
        catch ME
            throw(CORAerror("CORA:noops",SpS,obj));
        end
        res = aux_containsSpS_OuterSampling(obj, SpS, tol, maxEval);
        return;
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsPoint(SpS, p, tol)
    N = size(p,2);
    res = false([1 N]);
    
    persistent options
    if isempty(options)
        options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    end
    
    for i_N = 1:N
        G = SpS.G;
        c = SpS.c;
        [A0, Ai] = priv_getCoeffMatrices(SpS);
        
        m = size(G,2);
        
        beta = sdpvar(m,1,'full');
        s = sdpvar(size(A0,1),1,'full');
        t = sdpvar(1,1,'full');
        
        
        %  build constraint on matrix polynomial
        LMI = A0;
        
        % sum of matrix coefficients
        for i = 1:m
            LMI = LMI + beta(i) * Ai{i};
        end
        
        constraints = [...
            LMI + diag(s)>= 0,...
            G*beta + c == p(:,i_N),
            ];
        for i=1:size(A0,1)
            constraints = [constraints s(i)<=t -s(i)<=t];
        end
        
        cost = t;
        
        
        yalmipOptimizer = optimizer(constraints,cost,options,[],{beta,s,t});
        
        [sol,exitflag] = yalmipOptimizer();
        
        if exitflag == 1
            % The problem cannot be infeasible inherently, therefore
            % receiving exitflag = 1 means Sedumi encountered an internal
            % error and the point is within the set.
            sol{3} = 0;
        end

        res(i_N) = (sol{3}<=tol);
    end

end

function res = aux_containsSpS_InnerSampling(SpS1, SpS2, tol, N)
    
    if ~(dim(SpS1) == dim(SpS2))
        res = false;
        return;
    end

    persistent options
    if isempty(options)
        options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    end
    
    G = SpS1.G;
    m = size(G, 2);
    p = randPoint_(SpS1, N, 'uniform:billiardWalk');
    [A0,Ai] = priv_getCoeffMatrices(SpS1);
    
    for i=1:N
        % max lambda
        lambda = sdpvar(1, 1);
        
        % Randomize direction
        d = randn([dim(SpS1) 1]);
        d = d./norm(d);

        p0 = p(:,i);

        x = sdpvar(m, 1);

        constraints = [p0 + lambda*d == G*x];

        lmi = A0;
        for j=1:m
            lmi = lmi + Ai{j}*x(j);
        end

        constraints = [constraints, lmi>=0];
        
        yalmipOptimizer = optimizer(constraints, -lambda, options, [], {lambda, x});

        [sol, ~] = yalmipOptimizer();

        % identify a boundary point from p0
        lambda = sol{1};
        y = p0 + lambda * d;
    
        % if the boundary point is not in the outer SpS, the inner SpS
        % cannot be contained
        if ~SpS2.contains(y, 'exact', tol)
            res = false;
            return
        end
    end
    res = true;
end

function res = aux_containsSpS_OuterSampling(SpS1, SpS2, tol, N)
    if ~(dim(SpS1) == dim(SpS2))
        res = false;
        return;
    end

    persistent options
    if isempty(options)
        options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    end
    
    G = SpS2.G;

    c = center(SpS2);
    SpS1 = SpS1 - c;
    SpS2 = SpS2 - c;
    n = size(G, 1);
    m = size(G, 2);

    [A0,Ai] = priv_getCoeffMatrices(SpS2);

    p = randPoint_(SpS2, N, 'uniform:billiardWalk');

    for i=1:N
    % Compute random points on the boundary of SpS2
        % randomize direction
        d = randn([dim(SpS2) 1]);
        d = d./norm(d);

        p0 = p(:,i);

        % compute boundary point
        lambda = sdpvar(1,1);
        x = sdpvar(m, 1);

        constraints = [p0 + lambda * d == G*x];
        lmi = A0;
        for j=1:m
            lmi = lmi + Ai{j}*x(j);
        end
        constraints = [constraints, lmi>=0];

        yalmipOptimizer = optimizer(constraints, -lambda, options, [], {lambda, x});

        [sol, ~] = yalmipOptimizer();

        lambda = sol{1};

        y = p0 + lambda*d;

        % canonnical base vectors
        e=@(k,n) [zeros(k-1,1);1;zeros(n-k,1)];

        % compute the normal vector
        x = sdpvar(n, 1);
        Z = sdpvar(size(A0, 1), size(A0, 1));

        constraints = [trace(A0'*Z) <= 1];
        for j=1:m
            constraints = [constraints, trace(Ai{j}' * Z) == e(j,m)' * G' * x];
        end
        constraints = [constraints, Z>=0];
        objective = x' * y;

        yalmipOptimizer = optimizer(constraints,objective,options,[],{x, Z});

        [sol, ~] = yalmipOptimizer();

        eta = sol{1};
        
        % check the normal vector cutoff
        a = SpS1.supportFunc(eta./norm(eta));
        b = SpS2.supportFunc(eta./norm(eta));

        if a - tol > b
            res = false;
            return;
        end
    end
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
