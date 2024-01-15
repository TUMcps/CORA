function val = hausdorffDist(P,S)
% hausdorffDist - Calculates the Hausdorff distance between a polytope and
%    a set or a point; the cases are:
%    - identical sets (includes unbounded + unbounded, empty + empty): 0
%    - bounded + bounded (non-identical sets): bounded value > 0
%    - empty + bounded: Inf
%    - unbounded + bounded: Inf
%    - unbounded + unbounded (non-identical sets): Inf
%
% Syntax:
%    val = hausdorffDist(P,S)
%
% Inputs:
%    P - polytope object
%    S - contSet object of single point
%
% Outputs:
%    val - Hausdorff distance
%
% Examples:
%    A = [2 1; 0 2; -2 1; -1 -3; 1 -2];
%    b = ones(5,1);
%    P = polytope(A,b);
%    p = [4;4];
%   
%    val = hausdorffDist(P,p)
%
%    figure; hold on;
%    plot(P);
%    scatter(p(1),p(2));
%
% References: 
%    [1] S. Koenig, "Computational Aspects of the Hausdorff Distance in 
%        Unbounded Dimension", 2018
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/hausdorffDist

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       05-September-2018
% Last update:   16-December-2023 (MW, fix handling of equality constraints)
%                18-December-2023 (MW, fix special cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
val = 0;
% tolerance for property checks
tol = 1e-10;

% read out polytope object
[P,S] = findClassArg(P,S,'polytope');

persistent options
if isempty(options)
    if ~isSolverInstalled('mosek')
        options = optimoptions('quadprog','display','off');
    else
        options = [];
    end
end

% check dimensions
equalDimCheck(P,S);

% different cases for different types of sets
if isnumeric(S)

    % special cases (only linear programs) before regular computation
    % (facet and/or vertex enumeration and quadratic programs)
    if ~isBounded(P)
        val = Inf; return
    elseif representsa_(P,'emptySet',tol)
        if isempty(S)
            val = 0; return
        else
            val = Inf; return
        end
    end
    
    for i = 1:size(S,2)
        val_ = aux_distPolyPoint(P,S(:,i),options);
        val = max(val,val_);
    end
    
elseif isa(S,'polytope') || isa(S,'interval') || ...
       isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle')

    % special cases (only linear programs) before regular computation
    % (facet and/or vertex enumeration and quadratic programs)

    % empty cases
    if representsa_(P,'emptySet',tol)
        if representsa_(S,'emptySet',tol)
            val = 0; return
        else
            val = Inf; return
        end
    elseif representsa_(S,'emptySet',tol)
        val = Inf; return
    end
   
    % boundedness
    if ~isBounded(P)
        if isBounded(S)
            val = Inf; return
        elseif isequal(P,S,tol)
            val = 0; return
        else
            val = Inf; return
        end
    elseif ~isBounded(S)
        val = Inf; return
    end

    % convert set to polytope
    S = polytope(S);
    
    % compute distance d(P1,P2) = sup_{x \in P1} inf_{y \in P2) d(x,y)
    V = vertices(P);
    
    for i = 1:size(V,2)
       val_ = aux_distPolyPoint(S,V(:,i),options);
       val = max(val,val_);
    end
    
    % compute distance d(P2,P1) = sup_{x \in P2} inf_{y \in P1) d(x,y)
    V = vertices(S);
    
    for i = 1:size(V,2)
       val_ = aux_distPolyPoint(P,V(:,i),options);
       val = max(val,val_);
    end
   
else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',P,S));
end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_distPolyPoint(P,x,options)
% compute Hausdorff distance between a polytope and a single point
% according to Equation (7) in [1]

H = 2*eye(length(x));
f = -2*x;

[~,val] = quadprog(H,f,P.A,P.b,P.Ae,P.be,[],[],[],options);

val = sqrt(val + x'*x);

end

% ------------------------------ END OF CODE ------------------------------
