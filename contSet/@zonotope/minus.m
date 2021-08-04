function Z = minus(minuend,subtrahend,varargin)
% minus - overloaded '-' operator for approximating the Minkowski
%    difference of two zonotopes or a zonotope with a vector.
%    A - B = C <-> B + C \subseteq A
%
% Syntax:  
%    Z = minus(minuend,subtrahend)
%    Z = minus(minuend,subtrahend,type)
%    Z = minus(minuend,subtrahend,type,method)
%
% Inputs:
%    minuend - zonotope object
%    subtrahend - zonotope object or numerical vector
%    type - type of computation
%               - 'exact' (not implemented)
%               - 'approx' (default)
%    method - (optional) used algorithm
%               - 'althoff' (default)
%               - 'conZonotope'
%
% Outputs:
%    Z - Zonotpe after Minkowski difference
%
% Example: 
%    Z1 = zonotope([1 2 2; 0 0 2]);
%    Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
%
%    Z3 = Z1 - Z2;
%    Z4 = Z2 + Z3;
%
%    figure; hold on;
%    plot(Z1,[1 2], 'b');
%    plot(Z2,[1 2], 'r');
%    plot(Z3,[1 2], 'g');
%    plot(Z4,[1 2], 'k');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, conZonotope/minus

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      10-June-2015
% Last update:  22-July-2015
%               05-August-2015
%               20-August-2015
%               30-July-2016
%               14-November-2016
%               05-February-2021 (NK, added alternative algorithm)
%               06-May-2021 (MA, check added whether minuend is full dimensional)
% Last revision:---

%------------- BEGIN CODE --------------

% only center is shifted
if isnumeric(subtrahend)
    Z = minuend + (-subtrahend);
    return;
end

%check whether minuend is full dimensional
if isFullDim(minuend) 

    % enclose second set with zonotope
    subtrahend = zonotope(subtrahend);

    % parse input arguments
    method = 'althoff';

    if nargin > 2 && ~isempty(varargin{1})
        if strcmp(varargin{1},'exact')
           error(errNoExactAlg(cZ1,cZ2));
        end
    end
    if nargin > 3 && ~isempty(varargin{2})
       method = varargin{2}; 
    end

    % compute Minkowski difference with the approach from [1]
    if strcmp(method,'althoff')

        Z = minkDiffAlthoff(minuend,subtrahend);

    % compute Minkowski difference using constrained zonotopes    
    elseif strcmp(method,'conZonotope')

        Z = minkDiffConZono(minuend,subtrahend);

    else
        [id,msg] = errWrongInput('method');
        error(id,msg);
    end
    
else
    % display that minuend has to be full-dimensional
    disp('Minkowski difference only supports full-dimensional minuends')
    
    Z = [];
end

end


% Auxiliary Functions -----------------------------------------------------

function Z = minkDiffAlthoff(minuend,subtrahend)
% compute Minkowski difference using the approach in [1]

    %% determine generators to be kept
    % obtain halfspace representation
    [P,comb] = polytope(minuend);
    HorigTwice = get(P,'H');
    KorigTwice = get(P,'K');
    Horig = HorigTwice(1:0.5*end,:);
    
    % nr of subtrahend generators
    subtrahendGens = length(subtrahend.Z(1,:)) - 1;
    
    % intersect polytopes according to Theorem 3 of [1]
    delta_K = HorigTwice*subtrahend.Z(:,1);
    %delta_K = HorigTwice*(minuend.Z(:,1)-subtrahend.Z(:,1));
    for i = 1:subtrahendGens
        delta_K = delta_K + abs(HorigTwice*subtrahend.Z(:,i+1));
    end
    Korig_new = KorigTwice - delta_K;
    P_int = mptPolytope(HorigTwice,Korig_new);
    
    % remove redundant halfspaces and remember inices
    removedInd = removedHalfspaces(P_int,Horig);
    
    % Is the Minkowski difference an empty set?
    if removedInd == -inf
        % Minkowski difference is an empty set
        Z = [];
    else
        % Minkowski difference is not empty, but not all generators of the
        % minuend are required
        if ~isempty(removedInd)
            %% count generators that have been removed
            % nr of original generators
            gens = length(minuend.Z(1,2:end));
            % initialize indices 
            indices = zeros(1,gens);
            % add indices of generators that contributed to the removed
            % halfspace
            for i=1:length(removedInd)
                contributingIndices = comb(removedInd(i),:);
                indices(contributingIndices) = indices(contributingIndices) + 1;
            end
            % find generators to be removed
            n = length(minuend.Z(:,1));
            requiredRemovals = factorial(gens-1)/(factorial(n-2)*factorial(gens-n+1));  %binom(gens-1,dim-2);
            indRemove = (indices == requiredRemovals);

            %obtain reduced minuend
            G = minuend.Z(:,2:end);
            G(:,indRemove) = [];
            minuend.Z = [minuend.Z(:,1), G];

            %remove H, K
            C = Horig;
            C(removedInd,:) = [];
            d = Korig_new(1:0.5*end,:);
            d(removedInd) = [];
        else
            C = Horig;
            d = Korig_new(1:0.5*end,:);
        end
    
        %compute center
        c = minuend.Z(:,1) - subtrahend.Z(:,1);

        %obtain minuend generators
        G = minuend.Z(:,2:end);

        %% reverse computation from halfspace generation
        delta_d = d - C*minuend.Z(:,1) + C*subtrahend.Z(:,1);
        A_abs = abs(C*G);
        %alpha = A_abs\delta_d; %solve linear set of equations
        alpha = pinv(A_abs)*delta_d; %solve linear set of equations
        for i=1:length(alpha)
            Gnew(:,i) = alpha(i)*minuend.Z(:,i+1);
        end

        % instantiate Z
        Z = zonotope([c,Gnew]);
    end

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

    % compute zonotope under-approximation of the constrained zonotope
    Z = innerApprox(cZ,center(Z1));
end

function Z = innerApprox(cZ,c)
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
    
    Aeq = [cZ.Z(:,2:end)*T, zeros(dim(cZ),m_)];
    beq = c - cZ.Z(:,1) -cZ.Z(:,2:end)*p_;
    
    % construct objective function of the linear program
    f = -[zeros(1,m_),sum((cZ.Z(:,2:end)*T).^2,1)];
    
    % solve linear program to get interval inner-approximation of polytope
    options = optimoptions('linprog','display','off');
    
    x = linprog(f,A,b,Aeq,beq,[],[],options);
    
    c = x(1:m_); r = x(m_+1:end);
    int = interval(c-r,c+r);
    
    % compute transformation matrices
    off = p_ + T*center(int);
    S = T*diag(rad(int));
    
    % construct final zonotope
    c = cZ.Z(:,1) + cZ.Z(:,2:end)*off;
    G = cZ.Z(:,2:end)*S;
    
    Z = zonotope([c,G]);
end

%------------- END OF CODE --------------