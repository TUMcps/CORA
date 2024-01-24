function Z = minkDiff(minuend, subtrahend, varargin)
% minkDiff - computes the Minkowski difference of two zonotopes.
%       A - B = C <-> B + C \subseteq A
%
% Syntax:
%    Z = minkDiff(minuend,subtrahend)
%    Z = minkDiff(minuend,subtrahend,method)
%
% Inputs:
%    minuend - zonotope object
%    subtrahend - zonotope object or numerical vector
%    method - (optional) used algorithm
%     - 'approx' (default)
%     - 'exact' (only for 2d or aligned)
%     - 'inner' 
%     - 'outer'
%     - 'outer:coarse'
%     - 'outer:scaling' (subtrahend must be interval)
%     - 'inner:conZonotope'
%     - 'inner:RaghuramanKoeln' (implements [2])
%
% Outputs:
%    Z - zonotope after Minkowski difference
%
% Example:
%    Z1 = zonotope([1 2 2; 0 0 2]);
%    Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
%
%    Z3 = minkDiff(Z1,Z2);
%    Z4 = Z2 + Z3;
%
%    figure; hold on;
%    plot(Z1,[1 2], 'b');
%    plot(Z2,[1 2], 'r');
%    plot(Z3,[1 2], 'g');
%    plot(Z4,[1 2], 'k');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes",
%        arXiv, 2015.
%    [2] V. Raghuraman and J. P. Koeln. Set operations and order reductions
%        for constrained zonotopes. Automatica, 139, 2022. article no. 110204.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, conZonotope/minkDiff

% Authors:       Matthias Althoff, Niklas Kochdumper, Tobias Ladner
% Written:       10-June-2015
% Last update:   22-July-2015
%                05-August-2015
%                20-August-2015
%                30-July-2016
%                14-November-2016
%                05-February-2021 (NK, added alternative algorithm)
%                06-May-2021 (MA, check added whether minuend is full dimensional)
%                17-June-2022 (MA, case added that minuend and subtrahend are degenerate)
%                24-June-2022 (MA, under-approximation and over-approximation added)
%                01-July-2022 (MA, coarser and faster method for over-approximation added)
%                15-July-2022 (MA, method from Raghuraman & Koeln added)
%                27-July-2022 (MA, method from Raghuraman & Koeln rewritten for linprog)
%                09-November-2022 (MW, rename 'minkDiff', rename methods)
%                25-May-2023 (TL, added 'exact' for aligned zonotopes)
% Last revision: 25-May-2023 (TL, restructuring; more descriptive error messages)

% ------------------------------ BEGIN CODE -------------------------------

% list implemented algorithms
implementedAlgs = {'exact','inner','outer','outer:coarse','outer:scaling', ...
    'approx','inner:conZonotope','inner:RaghuramanKoeln'};

% parse inputs
if nargin < 2
    throw(CORAerror("CORA:notEnoughInputArgs",2))
elseif nargin > 3
    throw(CORAerror("CORA:tooManyInputArgs",3))
end

% check inputs
method = setDefaultValues({'approx'},varargin);
inputArgsCheck({ ...
    {minuend,'att','zonotope'}; ...
    {subtrahend,'att',{'contSet','numeric'}}; ...
    {method,'str',implementedAlgs}; ...
})

% check if subtrahend is numeric
if isnumeric(subtrahend)
    Z = minuend - subtrahend;
    return;
end

% check if dimensions match
n = dim(minuend);
if n ~= dim(subtrahend)
    throw(CORAerror('CORA:dimensionMismatch',minuend,subtrahend))
end

% check if subtrahend is zonotope
if ~isa(subtrahend,'zonotope')
    if ~(strcmp(method, 'outer:scaling') && isa(subtrahend,'interval'))
        warning(['CORA warning: zonotope/minkDiff: '...
            'Subtrahend is not a zonotope. Enclosing it with a zonotope.'])
    end
    % enclose second set with zonotope
    subtrahend = zonotope(subtrahend);
end

% check if subtrahend is point
if isempty(generators(subtrahend))
    Z = minuend - center(subtrahend);
    return
end

%check whether minuend is full dimensional
if isFullDim(minuend)

    % solution is exact for n==2 and enforced for this dimension [1,Prop.6]
    if n == 2
        method = 'exact';
    end

    % compute Minkowski difference with the approach from [1]
    if strcmpi(method, 'exact')
        if aux_areAligned(minuend, subtrahend)
            % exact solution for aligned sets according to [1,Prop.5]
            Z = zonotope(minuend.c - subtrahend.c, minuend.G - subtrahend.G);

        elseif n == 2
            % same method as 'inner' [1,Prop.6]
            Z = aux_minkDiffZono(minuend, subtrahend, method);

        else
            throw(CORAerror('CORA:wrongValue', 'third',...
                'No exact algorithm found: Sets have to be 2-dimensional or aligned.'))
        end

    elseif strcmp(method, 'outer') || strcmp(method, 'outer:coarse') || ...
            strcmp(method, 'inner') || strcmp(method, 'approx')
        Z = aux_minkDiffZono(minuend, subtrahend, method);

    elseif strcmp(method, 'inner:conZonotope')
        % compute Minkowski difference using constrained zonotopes
        Z = aux_minkDiffConZono(minuend, subtrahend);

    elseif strcmp(method, 'inner:RaghuramanKoeln')
        % compute Minkowski difference using [2]
        Z = aux_RaghuramanKoeln(minuend, subtrahend);

    elseif strcmp(method, 'outer:scaling')
        % compute  Minkowski difference using scaling
        Z = aux_minkDiffOuterInterval(minuend,subtrahend);
    end

else
    if isFullDim(subtrahend)
        % display that minuend has to be full-dimensional if subtrahend is
        % full-dimensional
        throw(CORAerror('CORA:wrongValue', 'first', 'full-dimensional zonotope'));
    else
        if rank(minuend) == rank(subtrahend)
            % transform the minuend and subtrahend into a space where the
            % minuend is full-dimensional using the singular value decomposition

            % range of minuend
            [U, S] = svd(generators(minuend));
            newDim = nnz(~all(withinTol(S,0))); % nr. of new dimensions
            P_minuend = U(1:newDim, :); % projection matrix
            
            % range of subtrahend
            [U, S] = svd(generators(subtrahend));
            newDim = nnz(~all(withinTol(S,0))); % nr. of new dimensions
            P_subtrahend = U(1:newDim, :); % projection matrix
            
            % Is the range of the minuend and subtrahend equal?
            if all(size(P_minuend) == size(P_subtrahend)) ...
                    && norm(P_minuend-P_subtrahend) <= 1e-10
                % project
                minuend_proj = P_minuend * minuend; % transformed minuend
                subtrahend_proj = P_minuend * subtrahend; % transformed subtrahend
                
                % solve problem in the transformed domain
                Z_proj = minkDiff(minuend_proj, subtrahend_proj, method);
                
                % transform solution back into the original domain
                Z = pinv(P_minuend) * Z_proj;

            else
                % no solution exists
                throw(CORAerror('CORA:wrongValue', 'first/second',...
                    ['for non full-dimensional zonotopes: ', ...
                    'projection matrix found by svd has to be equal']));
            end
        else
            % no solution exists
            throw(CORAerror('CORA:wrongValue', 'first/second',...
                ['for non full-dimensional zonotopes: ', ...
                'rank of generator matrix must be equal']));
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_minkDiffZono(minuend, subtrahend, method)
% compute Minkowski difference using the approach in [1]

%% determine generators to be kept
% obtain halfspace representation
P = polytope(minuend);
HorigTwice = P.A;
KorigTwice = P.b;
Horig = HorigTwice(1:0.5*end, :);

% nr of subtrahend generators
subtrahendGens = size(subtrahend.G,2);

% intersect polytopes according to Theorem 3 of [1]
delta_K = HorigTwice * subtrahend.c;
for i = 1:subtrahendGens
    delta_K = delta_K + abs(HorigTwice*subtrahend.G(:, i));
end
Korig_new = KorigTwice - delta_K;

C = Horig;
d = Korig_new(1:0.5*end, :);

%compute center
c = minuend.c - subtrahend.c;

%obtain minuend generators
G = minuend.G;

%% reverse computation from halfspace generation
n = dim(minuend);
if strcmp(method, 'inner') || (strcmp(method, 'exact') && n == 2)
    delta_d = d - C * minuend.c + C * subtrahend.c;
    A_abs = abs(C*G);
    dims = length(A_abs(1, :));
    % vector of cost function
    f = vecnorm(minuend.G, 2, 1);
    % A_abs x <= delta_d && x >= 0
    [alpha, ~, exitflag] = linprog(-f, [A_abs; -eye(dims)], [delta_d; zeros(dims, 1)]);
    if isempty(alpha) || exitflag ~= 1
        % return empty set with correct dimensions?
        Z = zonotope(zeros(n, 0));
        return
    end
    
elseif strcmp(method, 'outer') || strcmp(method, 'outer:coarse')
    % reduce delta_d using linear programming
    if strcmp(method, 'outer')
        d_shortened = aux_tightenHalfspaces(HorigTwice, Korig_new);
    else
        d_shortened = Korig_new;
    end

    % is set empty?
    if isempty(d_shortened)
        % return empty set with correct dimensions
        Z = zonotope(zeros(n, 0));
        return
    else
        % vector of cost function
        f = vecnorm(minuend.G, 2, 1);
        % obtain unrestricted A_abs and delta_d
        C = Horig;
        d = d_shortened(1:0.5*end, :);
        delta_d = d - C * minuend.c + C * subtrahend.c;
        A_abs = abs(C*G);
        dims = length(A_abs(1, :));
        % A_abs x >= delta_d && x >= 0
        [alpha, ~, exitflag] = linprog(f, [-A_abs; -eye(dims)], [-delta_d; zeros(dims, 1)]);
    end

elseif strcmp(method, 'approx')
    delta_d = d - C * minuend.c + C * subtrahend.c;
    A_abs = abs(C*G);
    % use pseudoinverse to compute an approximation
    alpha = pinv(A_abs) * delta_d; %solve linear set of equations using the pseudoinverse

else
    % should already be caught before
    throw(CORAerror('CORA:specialError',...
        sprintf("Unknown method: '%s'",method)))
end

% instantiate Z
Gnew = generators(minuend) * diag(alpha);
% remove all zero columns
Gnew = Gnew(:,~all(Gnew == 0,1));
Z = zonotope(c, Gnew);
end

function Z = aux_minkDiffOuterInterval(minuend,subtrahend)
% compute  Minkowski difference using scaling
% subtrahend must be an interval

if ~representsa_(subtrahend,'interval',eps)
    throw(CORAerror('CORA:wrongValue','second',...
        sprintf("interval (using method='outer:scaling')")))
end

% scale using interval enclosure
radMin = rad(interval(minuend));
radSub = rad(interval(subtrahend));
scale = 1 - radSub./radMin; % relative
Z = enlarge(minuend, scale); % outer

end

function Z = aux_minkDiffConZono(Z1, Z2)
% compute Minkowski difference based on constrained zonotopes

% convert first zonotope to constrained zonotope
cZ = conZonotope(Z1);

% compute Minkowski difference according to Theorem 1 in [1]
c = center(Z2);
G = generators(Z2);

cZ = cZ + (-c);

for i = 1:size(G, 2)
    cZ = and_(cZ+G(:, i), cZ+(-G(:, i)), 'exact');
end

% compute zonotope inner-approximation of the constrained zonotope
Z = aux_innerApprox(cZ);
end

function Z = aux_innerApprox(cZ)
% inner-approximate a constrained zonotope with a zonotope

% compute point satisfying all constraints with pseudo inverse
p_ = pinv(cZ.A) * cZ.b;

% compute null-space of constraints
T = null(cZ.A);

% transform boundary constraints of the factor hypercube
m = size(cZ.A, 2);
m_ = size(T, 2);

A = [eye(m); -eye(m)];
b = ones(2*m, 1);

A_ = A * T;
b_ = b - A * p_;

% construct constraint matrices for linear program
A = [A_, abs(A_*eye(m_))];
A = [A; zeros(m_), -eye(m_)];
b = [b_; zeros(m_, 1)];

% construct objective function of the linear program
f = -[zeros(1, m_), sum((cZ.G * T).^2, 1)];

% solve linear program to get interval inner-approximation of polytope
persistent options
if isempty(options)
    options = optimoptions('linprog', 'display', 'off');
end

[x, ~, exitflag] = linprog(f, A, b, [], [], [], [], options);

% check if constrained zonotope is empty
if isempty(x) || exitflag ~= 1
    % return empty set with correct dimensions
    Z = zonotope(zeros(dim(cZ), 0));
    return
end

c = x(1:m_);
r = x(m_+1:end);
r(r < 0) = 0;
int = interval(c-r, c+r);

% compute transformation matrices
off = p_ + T * center(int);
S = T * diag(rad(int));

% construct final zonotope
c = cZ.c + cZ.G * off;
G = cZ.G * S;

Z = zonotope(c, G);
end

function d_new = aux_tightenHalfspaces(C, delta_d)
% tighten halfspaces so that the polytope is identical with the same number
% of halfspaces

% loop over halfspaces
for i = 1:length(delta_d)
    % normal vector
    n = -C(i, :)';
    [~, d_new(i, 1), exitflag] = linprog(n, C, delta_d);
end
if exitflag ~= 1
    % linear program is infeasible since polytope is empty
    d_new = [];
else
    % values have the opposite sign
    d_new = -d_new;
end
end

function Z = aux_RaghuramanKoeln(Z_m, Z_s)
% Solves the Minkowski difference using the method described in
% [2, Theorem 7]. A direct implementation of [2, Theorem 7] can be
% found in the corresponding unit test. Here, we transform the linear
% program into the form required for linprog. A detailed derivation of the
% transformation can be found in the complementary documentation of CORA.

% Implemented by Matthias Althoff

% Computes the Minkowski difference according to [2]
% Z_m: minuend
% Z_s: subtrahend

% extract data
c_m = Z_m.center;
c_s = Z_s.center;
G_m = Z_m.generators;
G_s = Z_s.generators;
M = [G_m, G_s];

% dimension and nr of generators
n = length(c_m);
n_m = size(G_m, 2); % number of generators of Z_m
n_s = size(G_s, 2); % number of generators of Z_s

% create M_tilde
M_tilde = zeros(n*(n_m + n_s), n_m+n_s);
for i = 1:length(M(1, :))
    M_tilde(n*(i - 1)+1:n*(i - 1)+n, i) = M(:, i);
end


% linprog solved linear programs in the form (partially using LaTeX notation):
% \min_x f^T x
% such that:
% Ax <= b \\
% A_eq x = b_eq
% x_l <= x <= x_u

% A
a = kron(ones(1, n_m+2*n_s+1), eye(n_m));
I = eye(n_m*(n_m + 2 * n_s + 1));
A = [zeros(n_m, n_m+n_s), zeros(n_m, n_m*(n_m + 2 * n_s + 1)), a, zeros(n_m, n); ...
    zeros(n_m*(n_m + 2 * n_s + 1), n_m+n_s), I, -I, zeros(n_m*(n_m + 2 * n_s + 1), n); ...
    zeros(n_m*(n_m + 2 * n_s + 1), n_m+n_s), -I, -I, zeros(n_m*(n_m + 2 * n_s + 1), n)];

% b
b = [ones(n_m, 1); ...
    zeros(2*n_m*(n_m + 2 * n_s + 1), 1)];

% A_eq
A_eq = [M_tilde, -kron(eye(n_m+n_s), G_m), zeros(n*(n_m + n_s), n_m*n_s), zeros(n*(n_m + n_s), n_m), zeros(n*(n_m + n_s), n_m*(n_m + 2 * n_s + 1)), zeros(n*(n_m + n_s), n); ...
    zeros(n*n_s, n_m+n_s), zeros(n*n_s, n_m*(n_m + n_s)), -kron(eye(n_s), G_m), zeros(n*n_s, n_m), zeros(n*n_s, n_m*(n_m + 2 * n_s + 1)), zeros(n*n_s, n); ...
    zeros(n, n_m+n_s), zeros(n, n_m*(n_m + n_s)), zeros(n, n_m*n_s), -G_m, zeros(n, n_m*(n_m + 2 * n_s + 1)), -eye(n)];

% b_eq
b_eq = [zeros(n*(n_m + n_s), 1); ...
    reshape(-G_s, [], 1); ...
    c_s - c_m];

% f minimizes phi
f = [-ones(n_m+n_s, 1); zeros(2*n_m*(n_m + 2 * n_s + 1)+n, 1)];

% solve linear programming problem
[x, ~, exitflag] = linprog(f, A, b, A_eq, b_eq);

if exitflag ~= 1
    % no solution exists
    throw(CORAerror("CORA:specialError", 'No solution exists.'))
else
    % extract phi
    phi = x(1:n_m+n_s);
    % extract c_d
    c_d = x(end-n+1:end);
    % result
    Z = zonotope([c_d, [G_m, G_s] * diag(phi)]);
end

end

function res = aux_areAligned(minuend, subtrahend)
% check if generators are aligned

% extract generators
Gmin = generators(minuend);
Gsub = generators(subtrahend);

% check dimensions
if all(size(Gmin) == size(Gsub))
    % normalize
    normMin = max(vecnorm(Gmin));
    normSub = max(vecnorm(Gsub));

    if withinTol(normMin,0) && withinTol(normSub,0)
       % both have only all-zero generators
       res = true;
    elseif withinTol(normMin,0) || withinTol(normSub,0)
        % only one generator matrix has only all-zero generators
        res = false;
    else
       % normalize generators
       Gmin = Gmin / normMin;
       Gsub = Gsub / normSub;

       % generators have to be ordered
       res = compareMatrices(Gmin,Gsub,eps,'equal',true);
    end
else
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
