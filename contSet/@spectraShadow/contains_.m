function [res,cert,scaling] = contains_(SpS,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a spectrahedral shadow contains a set or a
%    point
%
% Syntax:
%    [res,cert,scaling] = contains_(SpS,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    SpS - spectraShadow object
%    obj - contSet object or vector
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       SpS will be detected as lying in SpS, which can be useful to
%       counteract errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in SpS, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in SpS).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(SpS - SpS.center) + SpS.center contains S.
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

    % Deal with trivial cases first
    if representsa(S, 'emptySet', tol)
        % Empty set is always contained
        res = true;
        cert = true;
        scaling = 0;
        return
    end

    % empty set
    if representsa(SpS, 'emptySet')
        if isnumeric(S)
            if isempty(S)
                res = true;
                scaling = 0;
                cert = true;
            else
                res = false;
                scaling = Inf;
                cert = true;
            end
        else
            if representsa(S, 'emptySet')
                res = true;
                scaling = 0;
                cert = true;
            else
                res = false;
                scaling = Inf;
                cert = true;
            end
        end
        return

        % fullspace
    elseif representsa(S, 'fullspace')
        if representsa(SpS, 'fullspace')
            res = true;
            cert = true;
            scaling = 0;
        else
            res = false;
            cert = true;
            scaling = inf;
        end
        return
    end

    % The code is not yet ready to deal with scaling
    cert = false;
    scaling = Inf;
    if scalingToggle
        throw(CORAerror('CORA:notSupported',...
                    "The computation of the scaling factor or cert " + ...
                    "for polynomial zonotopes is not yet implemented."));
    end

    if isnumeric(S)
        res = aux_containsPoint(SpS, S, tol);
        cert = true([1 size(S,2)]);
        return
    end

    switch method
        case 'exact'
            % Only works for polyhedral sets
            if isa(S, 'conZonotope') || isa(S, 'interval') || ...
               isa(S, 'polytope') || isa(S, 'zonoBundle') || ...
               isa(S, 'zonotope')

                V = vertices(S);
                [res, cert, scaling] = contains_(SpS, V, 'exact', tol, 0, certToggle, scalingToggle);
                res = all(res);
                cert = true;
                scaling = max(scaling);
            else
                throw(CORAerror("CORA:noExactAlg",SpS,S));
            end
        case 'approx'
            % Attempt to over-approximate S by an interval, and then use
            % the preceding method
            try
                S = interval(S);
            catch ME
                throw(CORAerror('CORA:noSpecificAlg',method,SpS,S));
            end
            V = vertices(S);
            [res, cert, scaling] = contains_(SpS, V, 'exact', tol, 0, certToggle, scalingToggle);
            res = all(res);
            cert = res;
            scaling = max(scaling);
        case 'sampling'
            % Try the sampling method, assuming we can somehow cast the set
            % as a spectrahedral shadow
            try
                S = spectraShadow(S);
            catch ME
                throw(CORAerror("CORA:noops",SpS,S));
            end
            res = aux_containsSpS_OuterSampling(S, SpS, tol, maxEval);
            cert = ~res;

        otherwise
            throw(CORAerror('CORA:noSpecificAlg',method,SpS,S));
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsPoint(SpS, p, tol)
    % checks if point is contained
    N = size(p,2);
    res = false([1 N]);
    
    % get options
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
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
        
        diagnostics = optimize(constraints,cost,options);

        t = value(t);
        
        if diagnostics.problem == 1
            % The problem cannot be infeasible inherently, therefore
            % receiving exitflag = 1 means Sedumi encountered an internal
            % error and the point is within the set.
            t = 0;
        end

        res(i_N) = (t<=tol);
    end

end

function res = aux_containsSpS_InnerSampling(SpS1, SpS2, tol, N)
        
    if ~(dim(SpS1) == dim(SpS2))
        res = false;
        return;
    end

    % get options
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
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
        
        diagnostics = optimize(constraints, -lambda, options);

        % identify a boundary point from p0
        lambda = value(lambda);
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

    % get options
    persistent options
    if isempty(options)
        if isSolverInstalled('mosek')
            options = sdpsettings('solver','mosek','verbose',0,'allownonconvex',0,'cachesolvers',1);
        else
            options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0,'cachesolvers',1);
        end
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

        diagnostics = optimize(constraints, -lambda, options);

        lambda = value(lambda);

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

        diagnostics = optimize(constraints,objective,options);

        eta = value(x);
        
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
