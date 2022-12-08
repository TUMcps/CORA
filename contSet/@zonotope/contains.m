function res = contains(Z,S,varargin)
% contains - determines if a zonotope contains a set or a point
%
% Syntax:  
%    res = contains(Z,S)
%    res = contains(Z,S,type)
%    res = contains(Z,S,type,tol)
%    res = contains(Z,S,type,tol,maxEval)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object or single point
%    type - (optional) method used for the containment check.
%       The available options are:
%           - 'exact' (default): Checks for containment by using either
%               'venum' or 'polymax', depending on the number of generators
%               of Z and the type of object S.
%           - 'approx': Checks for containment using 'st' (see below) if S 
%               is a zonotope, or any approximative method available
%               otherwise.
%           - 'venum': Checks for containment by enumerating all vertices
%               of S (see Algorithm 1 in [2]).
%           - 'polymax': Checks for containment by maximizing the
%               polyhedral norm w.r.t. Z over S (see Algorithm 2 in [2]).
%           - 'opt': Solves the containment problem via optimization
%               (see [2]) using the subroutine surrogateopt. If a solution
%               using 'opt' returns that Z1 is not contained in Z2, then
%               this is guaranteed to be the case. The runtime is
%               polynomial w.r.t. maxEval and the other inputs.
%           - 'st': Solves the containment problem using the
%               approximative method from [1]. If a solution using 'st'
%               returns that Z1 is contained in Z2, then this is guaranteed
%               to be the case. The runtime is polynomial w.r.t. all
%               inputs.
%    tol - (optional) tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors. Default is 100*eps.
%    maxEval - (optional) only, if 'opt' is used: Number of maximal
%       function evaluations for the surrogate optimization. By default,
%       this is set to max(500, 200 * number of generators of Z1).
%
% Outputs:
%    res - true/false
%
% Note: For low dimensions or number of generators, and if S is a point
% cloud with a very large number of points, it may be beneficial to first
% compute the halfspace representation of Z via Z.halfspace, and then call
% contains(Z,P).
%
% Example: 
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
% 
%    contains(Z1,Z2)
%    contains(Z1,Z3)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%    
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%    [2] A. Kulmburg, M. Althoff. "On the co-NP-Completeness of the
%        Zonotope Containment Problem", European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/contains, conZonotope/contains

% Author:       Matthias Althoff, Niklas Kochdumper, Adrian Kulmburg
% Written:      07-May-2007 
% Last update:  06-April-2017
%               14-September-2019
%               19-November-2019 (NK, changed to header format)
%               01-July-2021 (AK, modified input parsing, and implemented
%                               new methods from [2])
%               22-July-2022 (MA, method 'st' no longer requires YALMIP)
%               25-November-2022 (LS, method 'st' using sparse matrices)
%               25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------

    % pre-processing
    [resFound,vars] = pre_contains('zonotope',Z,S,varargin{:});
    
    % check premature exit
    if resFound
        % if result has been found, it is stored in the first entry of var
        res = vars{1}; return
    else
        Z = vars{1}; S = vars{2}; type = vars{3}; tol = vars{4}; maxEval = vars{5};
    end


    % init result
    res = true;
        
    % point (cloud) in zonotope containment (see [2], Definition 4)
    if isnumeric(S)
        % we return a logical array, one value per points
        res = false(1,size(S,2));

        % By comparing the complexity of computing the halfspace
        % representation and a linear program, one can see that computing
        % the halfspace representation is easier to compute when the
        % dimension n is <= 4, or when the number of generators m is <= n.
        % Thus, in those two cases, this is the easiest solution:
        if isempty(Z.halfspace) && (dim(Z) <= 4 || size(Z.generators,2) <= dim(Z))
            Z=halfspace(Z);
        end
        
        % If the halfspace-representation of the zonotope has already been
        % computed in advance, we can use this
        if ~isempty(Z.halfspace)
            for i = 1:size(S,2)
                %simple test: Is point inside the zonotope?
                inequality = Z.halfspace.H*S(:,i) < Z.halfspace.K | ...
                    withinTol(Z.halfspace.H*S(:,i),Z.halfspace.K,tol);
                res(i) = all(inequality);
            end
        else
            for i = 1:size(S,2)
                tmp = zonotopeNorm(Z,S(:,i)-Z.center);
                res(i) = tmp < 1 | withinTol(tmp,1,tol);
            end
        end
        
    % capsule/ellipsoid in zonotope containment
    elseif isa(S,'capsule') || isa(S,'ellipsoid')      
        P = mptPolytope(Z);
        res = contains(P,S);
        
    % taylm/polyZonotope in zonotope containment
    elseif isa(S,'taylm') || isa(S,'polyZonotope')
        if strcmp(type,'exact')
            throw(CORAerror('CORA:noops',S,Z));
        elseif strcmp(type,'approx')
            P = mptPolytope(Z);
            res = contains(P,S);
        else
            throw(CORAerror('CORA:noops',S,Z));
        end

    else
        % As of now, to respect backwards compatibility, the default method
        % for 'exact' is 'venum' if S is an interval, 'polymax'
        % otherwise. This will change in the future, once better heuristics
        % for these two methods are found, and their running time is then
        % analyzed precisely.
        
        % If S is an interval, we can transform it to a zonotope; if
        % 'exact' is chosen as a method, we change the method to 'venum'.
        if isa(S, 'interval')
            S = zonotope(S);
            if strcmp(type,'exact')
                type = 'venum';
            end
        end
        
        % Now, either S is a zonotope, and we can use the different
        % methods for zonotope containment, or it is something else. In the
        % latter case, we can send it to conZonotope.in, since it contains
        % all the instructions for the general case.
        if ~isa(S, 'zonotope')
            cZ = conZonotope(Z);
            % Here, we need to check that only the methods 'exact' or
            % 'approx' are used
            if ~strcmp(type,'exact') && ~strcmp(type,'approx')
                throw(CORAerror('CORA:notSupported',...
                    "The methods 'venum', 'polymax', 'st' and 'opt' are " + ...
                    "not yet implemented for containment checks other than zonotope-in-zonotope!"));
            end
            res = contains(cZ,S,type);
        else
            % We may now assume that S is a zonotope
            
            % Set adaptive maxEval, if needed
            if strcmp(type, 'opt') && maxEval == -1
                m1 = size(S.generators, 2);
                maxEval = max(500, 200 * m1);
            end
            
            % Simplify the zonotopes as much as possible, and rename them
            % to Z1 and Z2, to match the notation in [2]
            Z1 = deleteZeros(S);
            Z2 = deleteZeros(Z);
            
            % Depending on the method, choose the right subfunction
            switch type
                case {'exact', 'polymax'}
                    res = polyhedralMaximization(Z1, Z2, tol);
                case 'venum'
                    res = vertexEnumeration(Z1, Z2, tol);
                case {'approx', 'st'}
                    res = SadraddiniTedrake(Z1, Z2, tol);
                case 'opt'
                    % Check if surrogate optimization is available
                    GOT_installed = true;
                    try
                        optimoptions('surrogateopt');
                    catch
                        GOT_installed = false;
                    end
                    
                    if ~GOT_installed
                        [~, warn_id] = warning('You have not installed the Global Optimization Toolbox from MATLAB, and can therefore not use surrogateopt for solving the zonotope containment problem. Alternatively, the DIRECT algorithm will be used for now, but for improved results, please install the Global Optimization Toolbox.');
                        % Right after the first warning, disable it, so as
                        % to not clutter the command window
                        warning('off', warn_id);
                        res = DIRECTMaximization(Z1, Z2, tol, maxEval);
                    else
                        res = surrogateMaximization(Z1, Z2, tol, maxEval);
                    end
            end
        end
    end
end

function isIn = vertexEnumeration(Z1, Z2, tol)
% Solves the zonotope containment problem by checking whether the maximum
% value of the Z2-norm at one of the vertices of Z1 exceeds 1+tol. Checks
% every vertex until an answer is found (see also [2, Algorithm 1]).

% Let c1, c2 be the centers of Z1, Z2. We prepare the norm-function,
% returning the norm of v-c2 w.r.t. the Z2-norm, where v is a given vertex
% of Z1. Since v = c1 +- g_1 +- ... +- g_m, where g_i are the generators of
% Z1, the norm of v-c2 is the same as the norm of G*nu + c1-c2, where
% G is the generator matrix of Z1, nu = [+-1;...;+-1].
G = Z1.generators;
norm_Z2_nu = @(nu) zonotopeNorm(Z2, G*nu + Z1.center-Z2.center);

% Number of generators of Z1
m = size(G, 2);

% Create list of all combinations of generators we have to check (i.e., the
% choices of the +- signs from above). Note that this is a better strategy
% than computing the vertices directly, since it takes requires less
% memory.
% The next two lines produce all m-combinations of +-1
combinations = dec2bin(0:2^m-1)-'0';
combinations = 2*(combinations - 0.5);

for iter = combinations'
    if norm_Z2_nu(iter) > 1 + tol
        % If one vertex has norm larger than 1, it lies outside Z2, and
        % thus Z1 is not in Z2
        isIn = false;
        return
    end
end
% If no vertex of Z1 lies outside of Z2, Z1 lies in Z2
isIn = true;
end

function isIn = polyhedralMaximization(Z1, Z2, tol)
% Solves the zonotope containment problem by computing the maximal value of
% the polyhedral norm over Z1 w.r.t. Z2 (see also [2, Algorithm 2]).

% First, we need to shift Z2 so that it is centered around the origin
Z = Z2 - Z2.center;
% Then, we compute the halfspace representation, if it is not already
% computed
if isempty(Z.halfspace)
    Z = halfspace(Z);
end

% We then need the normalized halfspace-matrix
H_norm = Z.halfspace.H ./ Z.halfspace.K;

% Similarly to vertexEnum, we need to create combinations from the pool of
% vectors given by [G c1-c2], where G is the generator matrix of Z1, and
% c1, c2 are the centers of Z1, Z2.
V = [Z1.generators Z1.center-Z2.center];

% Polyhedral norm at a point p
poly_norm = @(p) max([0 max(H_norm * p)]);

% Store signs for the main step (these are the signs that matter for the
% decision of the sign of the x_j in [2, Algorithm 2]).
M = sign(H_norm * V);

n = size(M, 1);

% Iterate over each row of M
for i = 1:n
    mu = M(i,:); % Sign-combination
    maximum = poly_norm(V*mu'); % Compute the resulting polyhedral norm
    if maximum > 1+tol % If we found a point with polyhedral norm larger
                       % than 1, we can stop the algorithm.
        isIn = false;
        return
    end
end
% If no point has norm larger than 1, this means that Z1 is contained
% within Z2.
isIn = true;
end

function isIn = surrogateMaximization(Z1, Z2, tol, maxEval)
% Solves the zonotope containment problem by checking whether the maximum
% value of the Z2-norm at one of the vertices of Z1 exceeds 1+tol, using
% surrogate optimization, see also [2].

% Retrieve the generator matrix of Z1 and determine its size
G = Z1.generators;
m = size(G, 2);

% Prepare objective function to be minimized (see [2], or also
% vertexEnumeration above, since the idea is similar).
% Note that we need to use nu' here instead of nu, since for surrogateopt
% the points are given as 1xn-arrays, unlike fmincon for example.
norm_Z2_nu = @(nu) -zonotopeNorm(Z2, G*nu' + Z1.center-Z2.center);

% Setting up options
options = optimoptions('surrogateopt',...
    'Display', 'none',... % Suppress output
    'PlotFcn', [],... % Supress plot output
    'ObjectiveLimit', -1-tol,... % Stop when a maximum > 1+tol has been
                             ... % found (i.e., a point of Z1 outside
                             ... % of Z2 has been found)
    'MaxFunctionEvaluations', maxEval); % Set maximum number of
                                        % function evaluations
                                        
% Launch the optimization. Note that the fourth argument ensures that the
% points that are tested are integer-points, since a maximum can only
% happen on one of the vertices of Z1, meaning that nu should only have
% the values +-1, i.e., integer values.
[~, fval] = surrogateopt(norm_Z2_nu, -ones([m 1])', ones([m 1])', ones([m 1])', options);

if -fval > 1 + tol
    isIn = false;
else
    isIn = true;
end
    
end

function isIn = DIRECTMaximization(Z1, Z2, tol, maxEval)
% Solves the zonotope containment problem by checking whether the maximum
% value of the Z2-norm at one of the vertices of Z1 exceeds 1+tol, using
% DIRECT optimization, in case surrogate optimization is not available.

% Retrieve the generator matrix of Z1 and determine its size
G = Z1.generators;
m = size(G, 2);

% Prepare objective function to be minimized (see [2], or also
% vertexEnumeration above, since the idea is similar).
% Note that we need to use nu' here instead of nu, since for surrogateopt
% the points are given as 1xn-arrays, unlike fmincon for example.
Problem.f = @(nu) -zonotopeNorm(Z2, G*nu + Z1.center-Z2.center);

bounds = [-ones(m,1) ones(m,1)];
opts.maxevals = maxEval;
opts.showits = 0;
opts.limit = -1-tol;
opts.maxits    = inf;
                                        
% Launch the optimization. Note that the fourth argument ensures that the
% points that are tested are integer-points, since a maximum can only
% happen on one of the vertices of Z1, meaning that nu should only have
% the values +-1, i.e., integer values.

[fval,~,~] = Direct(Problem,bounds,opts);

if -fval > 1 + tol
    isIn = false;
else
    isIn = true;
end
    
end

function isIn = SadraddiniTedrake(Z1, Z2, tol)
% Solves the containment problem using the method described in
% [1, Corollary 4]. A direct implementation of [1, Corollary 4] can be 
% found in the corresponding unit test. Here, we transform the linear 
% progam into the form required for linprog. A detailed derivation of the 
% transformation can be found in the complementary documentation of CORA.   

% Implemented by Matthias Althoff, Lukas Schäfer

% extract data
c_s = Z1.center; % s: small, assuming that Z1 is contained in Z2
c_l = Z2.center; % l: large, assuming that Z1 is contained in Z2
G_s = Z1.generators; % s: small, assuming that Z1 is contained in Z2
G_l = Z2.generators; % l: large, assuming that Z1 is contained in Z2
n_s = size(G_s, 2); % nr of generators of small zonotope
n_l = size(G_l, 2); % nr of generators of large zonotope
n = size(G_s, 1); % system dimension

% linprog solved linear programs in the form (partially using LaTex notation):
% \min_x f^T x 
% such that:
% Ax <= b \\
% A_eq x = b_eq
% x_l <= x <= x_u 

% variable structure
idxVars.Gamma = 1:n_l*n_s; nVars = n_l*n_s;
idxVars.beta = nVars + (1:n_l); nVars = nVars + n_l;
idxVars.GammaAux = nVars + (1:n_l*n_s); nVars = nVars + n_l*n_s;
idxVars.betaAux = nVars + (1:n_l); nVars = nVars + n_l;

% equality constraint Aeq x = beq ----------------------------------------
Aeq.row = []; Aeq.col = []; Aeq.val = []; nCon = 0;
beq.row = []; beq.val = [];

% G_s = G_l*Gamma
tmp_Gamma = reshape(idxVars.Gamma,n_l,[]);
for ii = 1:n_s
    Aeq.row = [Aeq.row, repmat(nCon+(1:n),1,n_l)];
    Aeq.col = [Aeq.col, repelem(tmp_Gamma(:,ii)',1,n)];
    Aeq.val = [Aeq.val, G_l(:)'];
    beq.row = [beq.row, nCon+(1:n)];
    beq.val = [beq.val, G_s(:,ii)'];
    nCon = nCon + n;
end
% c_l - c_s = G_l*beta
Aeq.row = [Aeq.row, repmat(nCon+(1:n),1,n_l)];
Aeq.col = [Aeq.col, repelem(idxVars.beta,1,n)];
Aeq.val = [Aeq.val, G_l(:)'];
beq.row = [beq.row, nCon+(1:n)];
beq.val = [beq.val, (c_l-c_s)'];
Aeq.nCon = nCon + n;

% inequality constraint A x <= b -----------------------------------------
A.row = []; A.col = []; A.val = []; nCon = 0;
b.row = []; b.val = [];

% resolve absolute value: -Gamma - GammaAux <= 0
A.row = [A.row, repmat(nCon+(1:length(idxVars.Gamma)),1,2)];
A.col = [A.col, idxVars.Gamma,idxVars.GammaAux];
A.val = [A.val, ...
    -ones(1,length(idxVars.Gamma)),-ones(1,length(idxVars.GammaAux))];
nCon = nCon + length(idxVars.Gamma);
% resolve absolute value: Gamma - GammaAux <= 0
A.row = [A.row, repmat(nCon+(1:length(idxVars.Gamma)),1,2)];
A.col = [A.col, idxVars.Gamma,idxVars.GammaAux];
A.val = [A.val, ...
    ones(1,length(idxVars.Gamma)),-ones(1,length(idxVars.Gamma))];
nCon = nCon + length(idxVars.Gamma);
% resolve absolute value: -beta - betaAux <= 0
A.row = [A.row, repmat(nCon+(1:length(idxVars.beta)),1,2)];
A.col = [A.col, idxVars.beta,idxVars.betaAux];
A.val = [A.val, ...
    -ones(1,length(idxVars.beta)),-ones(1,length(idxVars.betaAux))];
nCon = nCon + length(idxVars.beta);
% resolve absolute value: beta - betaAux <= 0
A.row = [A.row, repmat(nCon+(1:length(idxVars.beta)),1,2)];
A.col = [A.col, idxVars.beta,idxVars.betaAux];
A.val = [A.val, ...
    ones(1,length(idxVars.beta)),-ones(1,length(idxVars.betaAux))];
nCon = nCon + length(idxVars.beta);

% sum([GammaAux,betaAux],2) <= 1+tol
tmp = [reshape(idxVars.GammaAux,n_l,[]),idxVars.betaAux'];
A.row = [A.row, repmat(nCon+(1:n_l),1,n_s+1)];
A.col = [A.col, tmp(:)'];
A.val = [A.val, ones(1,n_l*(n_s+1))];
b.row = [b.row, nCon+(1:n_l)];
b.val = [b.val, ones(1,n_l)*(1+tol)];
A.nCon = nCon + n_l;

% cost can be set arbitrarily --------------------------------------------
cost = ones(nVars,1);

% setup solver-specific problem data & solve -----------------------------

if isSolverInstalled('mosek')
    % merge constraint matrices
    prob.a = sparse([Aeq.row, Aeq.nCon+A.row],[Aeq.col, A.col], ...
        [Aeq.val, A.val],Aeq.nCon+A.nCon,nVars);
    % merge constraint boundaries
    prob.blc = [zeros(Aeq.nCon,1);-inf(A.nCon,1)];
    prob.blc(beq.row) = beq.val;
    prob.buc = zeros(Aeq.nCon+A.nCon,1);
    prob.buc(beq.row) = beq.val;
    prob.buc(Aeq.nCon+b.row) = b.val;
    % cost
    prob.c = cost';
    
    % solve linear programming problem
    [~,res] = mosekopt('minimize echo(0)',prob);
    if strcmp(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE') || ...
            strcmp(res.sol.bas.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
        isIn = true;
    else
        isIn = false;
    end

else
    % constraint matrices
    beq = sparse(beq.row,ones(1,length(beq.row)),beq.val,Aeq.nCon,1);
    Aeq = sparse(Aeq.row,Aeq.col,Aeq.val,Aeq.nCon,nVars);
    b = sparse(b.row,ones(1,length(b.row)),b.val,A.nCon,1);
    A = sparse(A.row,A.col,A.val,A.nCon,nVars);
    
    % solve linear programming problem
    [~,~,exitflag] = linprog(cost,A,b,Aeq,beq);
    
    isIn = (exitflag == 1);
end

end

%------------- END OF CODE --------------