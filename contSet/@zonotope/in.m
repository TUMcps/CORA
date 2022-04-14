function res = in(Z,obj,varargin)
% in - determines if obj is contained in Z
%
% Syntax:  
%    res = in(Z,obj)
%    res = in(Z,obj,type)
%    res = in(Z,obj,tol)
%    res = in(Z,obj,type,tol)
%    res = in(Z,obj,tol,type)
%    res = in(Z,obj,type,tol,maxEval)
%    res = in(Z,obj,tol,type,maxEval)
%
% Inputs:
%    Z - zonotope object
%    obj - contSet object or point
%    type - (optional) method used for the containment check.
%       The available options are:
%           - 'exact' (default): Checks for containment by using either
%               'venum' or 'polymax', depending on the number of generators
%               of Z and the type of object obj.
%           - 'approx': Checks for containment using 'st' if obj is a
%               zonotope, or any approximative method available otherwise.
%           - 'venum': Checks for containment by enumerating all vertices
%               of obj (see Algorithm 1 in [2]).
%           - 'polymax': Checks for containment by maximizing the
%               polyhedral norm w.r.t. Z over obj (see Algorithm 2 in [2]).
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
%    res - boolean whether obj is contained in Z, or not
%
% Note: For low dimensions or number of generators, and if obj is a point
% cloud with a very large number of points, it may be beneficial to first
% compute the halfspace representation of Z via Z.halfspace, and then call
% in(Z,P).
%
% Example: 
%    zono1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    zono2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    zono3 = zono2 + [3;0];
% 
%    in(zono1,zono2)
%    in(zono1,zono3)
% 
%    figure; hold on;
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'g');
%    
%    figure; hold on;
%    plot(zono1,[1,2],'b');
%    plot(zono3,[1,2],'r');
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
% See also: interval/in, conZonotope/in
%
% Author:       Matthias Althoff, Niklas Kochdumper, Adrian Kulmburg
% Written:      07-May-2007 
% Last update:  06-April-2017
%               14-Sep-2019
%               19-Nov-2019 (NK, changed to header format)
%               01-July-2021 (AK, modified input parsing, and implemented
%                               new methods from [2])
% Last revision:---

%------------- BEGIN CODE --------------

    res = true;

    % parse input arguments
    type = 'exact';
    tol = 100*eps;
    maxEval = -1;
    
    if nargin < 2
        error("Not enough input arguments.")
    elseif nargin >= 6
        error("Too many input arguments.")
    end
        
    if nargin >= 3
        if isstring(varargin{1}) || ischar(varargin{1})
            if ~ismember(varargin{1}, {'venum', 'polymax', 'exact', 'opt', 'approx', 'st'})
                error('Unknown method name.')
            else
                type = varargin{1};
            end
        elseif isscalar(varargin{1})
            if varargin{1} < 0
                error('tol must be non-negative.')
            else
                tol = varargin{1};
            end
        else
            error('Unknown parameter.')
        end
    end
    
    if nargin >= 4
        if isstring(varargin{2}) || ischar(varargin{2})
            if ~ismember(varargin{2}, {'venum', 'polymax', 'exact', 'opt', 'approx', 'st'})
                error('Unknown method name.')
            else
                type = varargin{2};
            end
        elseif isscalar(varargin{2})
            if varargin{2} < 0
                error('tol must be non-negative.')
            else
                tol = varargin{2};
            end
        else
            error('Unknown parameter.')
        end
    end
                
    if nargin >= 5
        maxEval = varargin{3};
        if ~isscalar(maxEval) || maxEval <= 0 || floor(maxEval) ~= maxEval
            error('steps must be a positive integer or inf.')
        end
    end

        
    % point (cloud) in zonotope containment (see [2], Definition 4)
    if isnumeric(obj)
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
            for i = 1:size(obj,2)
                %simple test: Is point inside the zonotope?
                N = length(Z.halfspace.K);
                inequality = (Z.halfspace.H*obj(:,i) - Z.halfspace.K <= tol * ones(N,1));
                res = (all(inequality));
            end
        else
            for i = 1:size(obj,2)
                res = zonotopeNorm(Z,obj(:,i)-Z.center) <= 1+tol;
                if res ~= 1
                    return; 
                end
            end
        end
        
    % capsule/ellipsoid in zonotope containment
    elseif isa(obj,'capsule') || isa(obj,'ellipsoid')      
        poly = mptPolytope(Z);
        res = in(poly,obj);
        
    % taylm/polyZonotope in zonotope containment
    elseif isa(obj,'taylm') || isa(obj,'polyZonotope')
        if strcmp(type,'exact')
            error('Exact containment check not yet possible for this set representation!');
        elseif strcmp(type,'approx')
            poly = mptPolytope(Z);
            res = in(poly,obj);
        else
            error("Containment check using a method other than 'exact' or 'approx' is not yet possible for this set representation!"); 
        end

    else
        % As of now, to respect backwards compatibility, the default method
        % for 'exact' is 'venum' if obj is an interval, 'polymax'
        % otherwise. This will change in the future, once better heuristics
        % for these two methods are found, and their running time is then
        % analyzed precisely.
        
        % If obj is an interval, we can transform it to a zonotope; if
        % 'exact' is chosen as a method, we change the method to 'venum'.
        if isa(obj, 'interval')
            obj = zonotope(obj);
            if strcmp(type,'exact')
                type = 'venum';
            end
        end
        
        % Now, either obj is a zonotope, and we can use the different
        % methods for zonotope containment, or it is something else. In the
        % latter case, we can send it to conZonotope.in, since it contains
        % all the instructions for the general case.
        if ~isa(obj, 'zonotope')
            cZ = conZonotope(Z);
            % Here, we need to check that only the methods 'exact' or
            % 'approx' are used
            if ~strcmp(type,'exact') && ~strcmp(type,'approx')
                error("The methods 'venum', 'polymax', 'st' and 'opt' are not yet implemented for containment checks other than zonotope-in-zonotope!")
            end
            res = in(cZ,obj,type);
        else
            % We may now assume that obj is a zonotope
            
            % Set adaptive maxEval, if needed
            if strcmp(type, 'opt') && maxEval == -1
                m1 = size(obj.generators, 2);
                maxEval = max(500, 200 * m1);
            end
            
            % Simplify the zonotopes as much as possible, and rename them
            % to Z1 and Z2, to match the notation in [2]
            Z1 = deleteZeros(obj);
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
% [1, Theorem 3]. Same notation as [1], except for Z1 and Z2 which are,
% respectively, the inbody and the circumbody.

% Implemented by Felix Gruber

% extract data
x = Z1.center;
y = Z2.center;
X = Z1.generators;
Y = Z2.generators;
nx = size(X, 2);
ny = size(Y, 2);

% yalmip constructs linear programming matrices
Gamma = sdpvar(ny, nx, 'full');
beda = sdpvar(ny, 1);
constraints = [...
    X == Y*Gamma, ...
    y - x == Y*beda, ...
    norm([Gamma, beda], Inf) <= 1+tol];
cost = []; % It suffices to check feasibility here, so no cost function
options = sdpsettings('solver','linprog', 'verbose',0, 'allownonconvex',0); 

% solve linear programming problem
yalmipOptimizer = optimizer(constraints, cost, options, [], {Gamma, beda});

[~, exitFlag] = yalmipOptimizer();

isIn = ~exitFlag;
end

%------------- END OF CODE --------------