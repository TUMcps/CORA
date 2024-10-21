function [P,varargout] = polytope(Z,varargin)
% polytope - converts a zonotope object to a polytope object
%
% Syntax:
%    P = polytope(Z)
%    P = polytope(Z,method)
%    [P,comb,isDeg] = polytope(Z,'exact')
%
% Inputs:
%    Z - zonotope object
%    method - approximation:
%               'exact': based on Theorem 7 of [1]
%               'outer:tight': uses interval outer-approximation
%               'outer:volume' uses volume
%
% Outputs:
%    P - polytope object
%    comb - (only method = 'exact') generator combinations corresponding to
%               the halfspaces
%    isDeg - (only method = 'exact') true/false whether polytope is
%               full-dimensional
%
% Example: 
%    Z = zonotope([1;-1],[3 2 -1; 1 -2 2]);
%    P = polytope(Z);
%
%    figure; hold on;
%    plot(P,[1,2],'r');
%    plot(Z,[1,2],'b');
%
% References:
%   [1] Althoff, M.; Stursberg, O. & Buss, M. Computing Reachable Sets
%       of Hybrid Systems Using a Combination of Zonotopes and Polytopes
%       Nonlinear Analysis: Hybrid Systems, 2010, 4, 233-249
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Matthias Althoff
% Written:       06-August-2018
% Last update:   10-November-2022 (MW, unify with other polytope functions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% fix a tolerance
tol = 1e-12;

% parse input arguments
method = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {method,'str',{'exact','outer:tight','outer:volume'}}});


if strcmp(method,'exact')
    % note: this method was previously called 'polytope'

    % obtain number of generators, dimensions
    Z = compact_(Z,'zeros',eps);
    comb = [];
    isDeg = ~isFullDim(Z,tol);

    if size(Z.G,2) == 0
        % generate equality constraint for the center vector
        n = dim(Z);
        C = [eye(n);-eye(n)]; d = [Z.c;-Z.c];        
    elseif ~isDeg
        [C,d,comb] = aux_polytope_fullDim(Z);
    else
        [C,d] = aux_polytope_degenerate(Z,tol);
    end

    % other return values
    varargout{1} = comb;
    varargout{2} = isDeg;

    % instantiate polytope
    P = polytope(C,d);

elseif startsWith(method,'outer')
    P = aux_polytope_outer(Z,method);

end

% polytope is definitely bounded
P.bounded.val = true;

end


% Auxiliary functions -----------------------------------------------------

function [C,d,comb] = aux_polytope_fullDim(Z)

c = Z.c; G = Z.G;
[n,nrGen] = size(G);

if n == 1
    C = 1;
    comb = [];
else
    % get number of possible facets
    comb = combinator(nrGen,n-1,'c');
    % bypass bug in combinator (rows with all-zeros?!)
    comb = comb(any(comb,2),:);
    nrComb = size(comb,1);

    % build C matrices for inequality constraint C*x < d
    C = zeros(nrComb,n);
    for i=1:nrComb
        % compute n-dimensional cross product with each combination
        C(i,:) = ndimCross(G(:,comb(i,:)));
    end
    % normalize each normal vector
    C = C ./ vecnorm(C',2,1)';

    % remove NaN rows due to rank deficiency
    index = find(sum(isnan(C),2));
    if ~isempty(index)
        C(index,:) = [];
    end
end

% determine offset vector in addition to center
deltaD = sum(abs(C*G),2);
 
% construct the overall inequality constraints
d = [C*c + deltaD; -C*c + deltaD];
C = [C; -C];

end

function [C,d] = aux_polytope_degenerate(Z,tol)

% read out information about zonotope
c = Z.c; G = Z.G;
[n,nrGen] = size(G);

% singular value decomposition
[U,S,~] = svd(G);
S = [S,zeros(n,n-nrGen)];

% state space transformation
Z_ = U'*[c,G];

% remove dimensions with all zeros
ind = find(diag(S) <= tol);
ind_ = setdiff(1:size(S,1),ind);

if ~isempty(ind)
    % compute polytope in transformed space
    P = polytope(zonotope(Z_(ind_,:)));

    % transform back to original space
    A = [P.A,zeros(size(P.A,1),length(ind))]*U';
    b = P.b;

    % add equality constraint restricting polytope to null-space
    C = [A;U(:,ind)';-U(:,ind)'];
    d = [b;U(:,ind)'*c;-U(:,ind)'*c];
end

end

function P = aux_polytope_outer(Z)
% note: this method was previously called 'enclosingPolytope'

if strcmp(method,'outer:tight')
    % solution1 (axis-aligned):
    Zred = zonotope(interval(Z));
    P = polytope(Zred); 
    % solution 2 (method C):
    % Zred = reduce(Z,'methC',1,filterLength);
    Zred = reduce(Z,'pca');
    Zred = aux_repair(Zred,Z);
    Padd = polytope(Zred);
    % intersect results
    P = P & Padd;
elseif strcmp(method,'outer:volume')
    % solution 1 (method C):
    % Zred1 = reduce(Z,'methC',1,filterLength);
    Zred1 = reduce(Z,'pca');
    Zred1 = aux_repair(Zred1,Z);
    vol1 = volume(Zred1);
    % solution2 (axis-aligned):
    Zred2 = zonotope(interval(Z));
    Zred2 = aux_repair(Zred2,Z);
    vol2 = volume(Zred2);

    if vol1 < vol2
        P = polytope(Zred1);
    else 
        P = polytope(Zred2);
    end
end

end

function Zrep = aux_repair(Z,Zorig)
% repair zonotope if there is no length in one dimenion

% get length of each dimension
len = 2*rad(interval(Z));

% find zero lengths
index = find(len==0);

if ~isempty(index)
    %construct zonotope to be added
    origLen = 2*rad(interval(Zorig));
    origCenter = center(interval(Zorig));
    
    %get Zmatrix
    Zmat = [Z.c,Z.G];
    
    for i=1:length(index)
        
        ind = index(i);
        %replace center
        Zmat(ind,1) = origCenter(ind);
        %replace generator value
        Zmat(ind,ind+1) = 0.5*origLen(ind);
    end

    % instantiate zonotopes
    Zrep = zonotope(Zmat);
else
    Zrep = Z;
end

end
    
% ------------------------------ END OF CODE ------------------------------
