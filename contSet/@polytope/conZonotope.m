function cZ = conZonotope(P,varargin)
% conZonotope - converts a polytope object into a constrained zonotope;
%    implementation according to Theorem 1 from [1]
%
% Syntax:
%    cZ = conZonotope(P)
%    cZ = conZonotope(P,method)
%    cZ = conZonotope(P,method,box)
%
% Inputs:
%    P - polytope object
%    method - (optional) conversion method
%             - 'exact:supportFunc': using support functions (default)
%             - 'exact:vertices': using vertex enumeration
%    box - (optional) box outer approximation of P
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    A = [-1 0;0 -1;1 1];
%    b = [1;1;1];
%    P = polytope(A,b);
%    cZ = conZonotope(P);
%
%    figure; hold on
%    plot(cZ);
%    plot(P,[1,2],'r--');
%
% References: 
%    [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%        estimation and fault detection"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       13-May-2018
% Last update:   28-April-2019 (MA, code shortened)
%                13-December-2022 (MW, add support-function based method)
%                14-July-2023 (MW, add support for empty polytopes)
%                27-July-2023 (MW, incorporate equality constraints)
%                03-January-2024 (MW, speed up supportFunc method, handle unbounded)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% dimension
n = dim(P);

% set default method
[method,B] = setDefaultValues({'exact:supportFunc',...
    interval(zeros(n,0),zeros(n,0))},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'}, ...
                {method,'str',{'exact:vertices','exact:supportFunc'}}, ...
                {B,'att','interval'}});

% number of constraints
nrIneq = size(P.A,1);
nrEq = size(P.Ae,1);

if representsa_(P,'fullspace',0)
    % conversion of fullspace object not possible
    throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
        'can therefore not be converted into a constrained zonotope.']));

elseif strcmp(method,'exact:vertices')

    % calculate the vertices of the polytope (check for unboundedness)
    try
        V = vertices(P);
    catch ME
        if ~isempty(P.bounded.val) && ~P.bounded.val
             throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
                'can therefore not be converted into a constrained zonotope.']));
        end
        rethrow(ME);
    end

    % no vertices -> empty set
    if isempty(V)
        cZ = conZonotope.empty(n);
        return
    end

    % read out all constraints
    A_all = [P.A; P.Ae];
    b_all = [P.b; P.be];
    
    % calculate a bounding box for the constrained zonotope
    minV = min(V,[],2);
    maxV = max(V,[],2);
    
    % compute center and generator matrix
    c = 0.5 * (maxV + minV);
    G = diag(0.5 * (maxV - minV));
    
    % Calculate the lower bound sigma for A*x \in [sigma,b] (Thm. 1 in [1])
    sigma = min(A_all*V,[],2);
    
    % Construct constrained zonotope object according to eq. (21) in [1]
    G_ = [G, zeros(size(G,1),nrIneq+nrEq)];
    A_ = [A_all*G, diag((sigma-b_all)./2)];
    b_ = (b_all+sigma)./2 - A_all*c;

elseif strcmp(method,'exact:supportFunc')
    
    % compute bounding box
    B = interval(P);

    % check if P is unbounded
    if ~P.bounded.val
        % conversion of unbounded object not possible
        throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
            'can therefore not be converted into a constrained zonotope.']));
    elseif P.emptySet.val
        cZ = conZonotope.empty(n);
        return
    end
    
    % compute center and generators
    c = center(B);
    G = diag(0.5 * (supremum(B) - infimum(B)));

    % read out constraints of polytope
    A_all = [P.A; P.Ae];
    b_all = [P.b; P.be];
    
    % compute lower bound in the direction of halfspaces
    sigma = zeros(nrIneq,1);
    for a=1:nrIneq
        sigma(a) = supportFunc_(B,A_all(a,:)','lower');
        % any lower bound is Inf -> polytope is empty
        if sigma(a) == Inf
            cZ = conZonotope.empty(n); return
        end
    end
    % same for equality constraints
    if ~isempty(P.Ae)
        % no need to compute the value
        sigma = [sigma; P.be];
    end
    
    % Construct constrained zonotope object according to eq. (21) in [1]
    G_ = [G, zeros(size(G,1),nrIneq+nrEq)];
    A_ = [A_all*G, diag((sigma-b_all)./2)];
    b_ = (b_all+sigma)./2 - A_all*c;

end

% init constained zonotope
cZ = conZonotope(c,G_,A_,b_);

% ------------------------------ END OF CODE ------------------------------
