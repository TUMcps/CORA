function val = hausdorffDist(P,S)
% hausdorffDist - Calculates the Hausdorff distance between a polytope and
%    a set or a point
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
% See also: quadZonotope/hausdorffDist

% Authors:       Niklas Kochdumper
% Written:       05-September-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init result
val = 0;

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

% differnt cases for different types of sets
if isnumeric(S)
    
    for i = 1:size(S,2)
        val_ = aux_distPolyPoint(P,S(:,i),options);
        val = max(val,val_);
    end
    
elseif isa(S,'polytope') || isa(S,'interval') || ...
       isa(S,'zonotope') ||  isa(S,'conZonotope') || ...
       isa(S,'zonoBundle')
   
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
    error(noops(P,S));
end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_distPolyPoint(P,x,options)
% compute Hausdorff distance between a polytope and a single point
% according to Equation (7) in [1]

H = 2*eye(length(x));
f = -2*x;

A = P.A; 
b = P.b;

[~,val] = quadprog(H,f,A,b,[],[],[],[],[],options);

val = sqrt(val + x'*x);

end

% ------------------------------ END OF CODE ------------------------------
