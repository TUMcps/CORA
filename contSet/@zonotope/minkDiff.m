function Z = minkDiff(minuend,subtrahend,varargin)
% minkDiff - computes he Minkowski difference of two zonotopes.
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
%               - 'inner' (default)
%               - 'outer' 
%               - 'outer:coarse'
%               - 'approx'
%               - 'inner:conZonotope'
%               - 'inner:RaghuramanKoeln' (implements [2])
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      10-June-2015
% Last update:  22-July-2015
%               05-August-2015
%               20-August-2015
%               30-July-2016
%               14-November-2016
%               05-February-2021 (NK, added alternative algorithm)
%               06-May-2021 (MA, check added whether minuend is full dimensional)
%               17-June-2022 (MA, case added that minuend and subtrahend are degenerate)
%               24-June-2022 (MA, under-approximation and over-approximation added)
%               01-July-2022 (MA, coarser and faster method for over-approximation added)
%               15-July-2022 (MA, method from Raghuraman & Koeln added)
%               27-July-2022 (MA, method from Raghuraman & Koeln rewritten for linprog)
%               09-November-2022 (MW, rename 'minkDiff', rename methods)
% Last revision:---

%------------- BEGIN CODE --------------

% only center is shifted
if isnumeric(subtrahend)
    Z = minuend + (-subtrahend);
    return;
end



% list implemented algorithms
implementedAlgs = {'inner','outer','outer:coarse','approx',...
    'inner:conZonotope','inner:RaghuramanKoeln'};

%check whether minuend is full dimensional
if isFullDim(minuend)

    % enclose second set with zonotope
    subtrahend = zonotope(subtrahend);

    %% parse input arguments
    % default value
    method = 'approx';
    
    % read optional argument
    if nargin > 2 && ~isempty(varargin{1})
        if ~any(strcmp(varargin{1},implementedAlgs))
            throw(CORAerror('CORA:wrongValue','third',...
                strjoin(implementedAlgs,', ')));
        end
        method = varargin{1};
    end

    % compute Minkowski difference with the approach from [1]
    if strcmp(method,'outer') || strcmp(method,'outer:coarse') || ...
            strcmp(method,'inner') || strcmp(method,'approx')

        Z = minkDiffZono(minuend,subtrahend,method);

    % compute Minkowski difference using constrained zonotopes    
    elseif strcmp(method,'inner:conZonotope')
        Z = minkDiffConZono(minuend,subtrahend);
    % compute Minkowski difference using [2]   
    elseif strcmp(method,'inner:RaghuramanKoeln')
        Z = RaghuramanKoeln(minuend, subtrahend);
    end
    
else
    if isFullDim(subtrahend)
        % display that minuend has to be full-dimensional if subtrahend is
        % full-dimensional
        throw(CORAerror('CORA:wrongValue','first','full-dimensional zonotope'));
    else
        if rank(minuend)==rank(subtrahend)
            % transform the minuend and subtrahend into a space where the
            % minuend is full-dimensional using the singular value decomposition
            % range of minuend
            [U,S] = svd(generators(minuend));
            newDim = length(S(1,:)); % new dimension
            P_minuend = U(1:newDim,:); % projection matrix
            % range of subtrahend
            [U,S] = svd(generators(subtrahend));
            newDim = length(S(1,:)); % new dimension
            P_subtrahend = U(1:newDim,:); % projection matrix
            % Is the range of the minuend and subtrahend equal?
            if norm(P_minuend - P_subtrahend) <= 1e-10
                minuend_proj = P_minuend*minuend; % transformed minuend
                subtrahend_proj = P_minuend*subtrahend; % transformed subtrahend
                % solve problem in the transformed domain
                Z_proj = minkDiff(minuend_proj,subtrahend_proj,varargin);
                % transform solution back into the original domain
                Z = P_minuend'*Z_proj;
            else
                % no solution exists --> replace by error 
                Z = [];
            end
        else
            % no solution exists --> replace by error 
            Z = [];
        end
    end
end

end


% Auxiliary Functions -----------------------------------------------------

function Z = minkDiffZono(minuend,subtrahend,method)
% compute Minkowski difference using the approach in [1]

    %% determine generators to be kept
    % obtain halfspace representation
    P = mptPolytope(minuend);
    HorigTwice = get(P,'H');
    KorigTwice = get(P,'K');
    Horig = HorigTwice(1:0.5*end,:);
    
    % nr of subtrahend generators
    subtrahendGens = length(subtrahend.Z(1,:)) - 1;
    
    % intersect polytopes according to Theorem 3 of [1]
    delta_K = HorigTwice*subtrahend.Z(:,1);
    for i = 1:subtrahendGens
        delta_K = delta_K + abs(HorigTwice*subtrahend.Z(:,i+1));
    end
    Korig_new = KorigTwice - delta_K;
    
    C = Horig;
    d = Korig_new(1:0.5*end,:);

    %compute center
    c = minuend.Z(:,1) - subtrahend.Z(:,1);

    %obtain minuend generators
    G = minuend.Z(:,2:end);

    %% reverse computation from halfspace generation
    % approx is exact for dim==2 and enforced for this dimension
    if strcmp(method,'inner') || (length(c) == 2)
        delta_d = d - C*minuend.Z(:,1) + C*subtrahend.Z(:,1);
        A_abs = abs(C*G);
        dim = length(A_abs(1,:));
        % vector of cost function
        f = vecnorm(minuend.Z(:,2:end),2,1);
        % A_abs x <= delta_d && x >= 0
        [alpha,~,exitflag] = linprog(-f,[A_abs; -eye(dim)],[delta_d;zeros(dim,1)]);
        if isempty(alpha) || exitflag ~= 1
            Z = [];
            return 
        end
    elseif strcmp(method,'outer') || strcmp(method,'outer:coarse') 
        % reduce delta_d using linear programming
        if strcmp(method,'outer')
            d_shortened = tightenHalfspaces(HorigTwice,Korig_new);
        else
            d_shortened = Korig_new;
        end
        % is set empty?
        if ~isempty(d_shortened)
            % vector of cost function
            f = vecnorm(minuend.Z(:,2:end),2,1);
            % obtain unrestricted A_abs and delta_d
            C = Horig;
            d = d_shortened(1:0.5*end,:);
            delta_d = d - C*minuend.Z(:,1) + C*subtrahend.Z(:,1);
            A_abs = abs(C*G);
            dim = length(A_abs(1,:));
            % A_abs x >= delta_d && x >= 0
            [alpha,~,exitflag] = linprog(f,[-A_abs; -eye(dim)],[-delta_d;zeros(dim,1)]);
        else
            Z = [];
            return 
        end
    elseif strcmp(method,'approx') 
        delta_d = d - C*minuend.Z(:,1) + C*subtrahend.Z(:,1);
        A_abs = abs(C*G);
        % use pseudoinverse to compute an approximation
        alpha = pinv(A_abs)*delta_d; %solve linear set of equations using the pseudoinverse
    end
    
    % instantiate Z
    Gnew = generators(minuend)*diag(alpha);
    Z = zonotope(c,Gnew);
end

function Z = minkDiffConZono(Z1,Z2)
% compute Minkowski difference based on constrained zonotopes

    % convert first zonotope to constrained zonotope
    cZ = conZonotope(Z1);
    
    % compute Minkowski difference according to Theorem 1 in [1]
    c = center(Z2);
    G = generators(Z2);
    
    cZ = cZ + (-c);
    
    for i = 1:size(G,2)
        cZ = (cZ + G(:,i)) & (cZ + (-G(:,i)));
    end

    % compute zonotope inner-approximation of the constrained zonotope
    Z = innerApprox(cZ);
end

function Z = innerApprox(cZ)
% inner-approximate a constrained zonotope with a zonotope

    % compute point satisfying all constraints with pseudo inverse
    p_ = pinv(cZ.A)*cZ.b;
    
    % compute null-space of constraints
    T = null(cZ.A);
    
    % transform boundary constraints of the factor hypercube
    m = size(cZ.A,2);
    m_ = size(T,2);
    
    A = [eye(m);-eye(m)];
    b = ones(2*m,1);
    
    A_ = A*T;
    b_ = b - A*p_;
    
    % construct constraint matrices for linear program
    A = [A_, abs(A_*eye(m_))];
    A = [A; zeros(m_) -eye(m_)];
    b = [b_; zeros(m_,1)];
    
    % construct objective function of the linear program
    f = -[zeros(1,m_),sum((cZ.Z(:,2:end)*T).^2,1)];
    
    % solve linear program to get interval inner-approximation of polytope
    options = optimoptions('linprog','display','off');
    
    [x,~,exitflag] = linprog(f,A,b,[],[],[],[],options);
    
    % check if constrained zonotope is empty
    if isempty(x) || exitflag ~= 1
        Z = [];
        return 
    end
        
    c = x(1:m_); r = x(m_+1:end); r(r<0) = 0;
    int = interval(c-r,c+r);
    
    % compute transformation matrices
    off = p_ + T*center(int);
    S = T*diag(rad(int));
    
    % construct final zonotope
    c = cZ.Z(:,1) + cZ.Z(:,2:end)*off;
    G = cZ.Z(:,2:end)*S;
    
    Z = zonotope(c,G);
end

function d_new = tightenHalfspaces(C,delta_d)
% tighten halfspaces so that the polytope is identical with the same number
% of halfspaces

    % loop over halfspaces 
    for i = 1:length(delta_d)
        % normal vector
        n = -C(i,:)';
        [~,d_new(i,1),exitflag] = linprog(n,C,delta_d);
    end
    if exitflag ~= 1
        % linear program is infeasible since polytope is empty
        d_new = []; 
    else
        % values have the opposite sign
        d_new = -d_new;
    end
end

function Z = RaghuramanKoeln(Z_m, Z_s)
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
M_tilde = zeros(n*(n_m + n_s), n_m + n_s);
for i=1:length(M(1,:))
    M_tilde(n*(i-1)+1:n*(i-1)+n,i) = M(:,i);
end


   
% linprog solved linear programs in the form (partially using LaTeX notation):
% \min_x f^T x 
% such that:
% Ax <= b \\
% A_eq x = b_eq
% x_l <= x <= x_u 

% A
a = kron(ones(1,n_m+2*n_s+1),eye(n_m));
I = eye(n_m*(n_m+2*n_s+1));
A = [zeros(n_m, n_m+n_s), zeros(n_m, n_m*(n_m+2*n_s+1)), a, zeros(n_m, n); ...
    zeros(n_m*(n_m+2*n_s+1), n_m+n_s), I, -I, zeros(n_m*(n_m+2*n_s+1),n); ...
    zeros(n_m*(n_m+2*n_s+1), n_m+n_s),-I, -I, zeros(n_m*(n_m+2*n_s+1),n)];

% b
b = [ones(n_m,1); ...
    zeros(2*n_m*(n_m+2*n_s+1),1)];

% A_eq
A_eq = [M_tilde, -kron(eye(n_m+n_s), G_m), zeros(n*(n_m+n_s),n_m*n_s), zeros(n*(n_m+n_s), n_m), zeros(n*(n_m+n_s), n_m*(n_m+2*n_s+1)), zeros(n*(n_m+n_s), n); ...
    zeros(n*n_s, n_m + n_s), zeros(n*n_s, n_m*(n_m+n_s)), -kron(eye(n_s),G_m), zeros(n*n_s,n_m), zeros(n*n_s, n_m*(n_m+2*n_s+1)), zeros(n*n_s,n); ...
    zeros(n, n_m+n_s), zeros(n, n_m*(n_m + n_s)), zeros(n, n_m*n_s), -G_m, zeros(n, n_m*(n_m+2*n_s+1)), - eye(n)];

% b_eq
b_eq = [zeros(n*(n_m+n_s),1); ...
    reshape(-G_s,[],1); ...
    c_s-c_m];

% f minimizes phi
f = [-ones(n_m+n_s,1); zeros(2*n_m*(n_m+2*n_s+1) + n,1)];

% solve linear programming problem
[x,~,exitflag] = linprog(f,A,b,A_eq,b_eq);

if exitflag ~= 1
    % no solution exists
    Z = [];
else
    % extract phi
    phi = x(1:n_m+n_s);
    % extract c_d
    c_d = x(end-n+1:end);
    % result
    Z = zonotope([c_d,[G_m, G_s]*diag(phi)]);
end

end

%------------- END OF CODE --------------