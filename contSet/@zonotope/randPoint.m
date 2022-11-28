function p = randPoint(Z,varargin)
% randPoint - generates random points within a zonotope
%
% Syntax:  
%    p = randPoint(Z)
%    p = randPoint(Z,N)
%    p = randPoint(Z,N,type)
%    p = randPoint(Z,'all','extreme')
%    p = randPoint(Z,N,'gaussian',pr)
%
% Inputs:
%    Z - zonotope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', or 'gaussian')
%    pr - probability that a value is within the set (only type = 'gaussian')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z = zonotope([1;0],rand(2,5));
%    p = randPoint(Z);
% 
%    plot(Z); hold on;
%    scatter(p(1,:),p(2,:),16,'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/randPoint

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       23-September-2008 
% Last update:   25-June-2021 (MP, add type gaussian)
%                19-August-2022 (MW, integrate standardized pre-processing)
% Last revision: ---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_randPoint('zonotope',Z,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of vars
    p = vars{1}; return
else
    % assign variables
    Z = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
    % sampling with 'gaussian' is done in contSet method
    if strcmp(type,'gaussian')
        p = randPoint@contSet(Z,N,type,pr); return
    end
end

% get object properties
c = center(Z); G = generators(Z); n = dim(Z);
        
% no generators: zonotope is just a point
if isempty(G) || ~any(any(G))
    % replicate center N times
    p = repmat(c,1,N); return
end
        
% generate different types of random points
if strcmp(type,'standard')
    % degeneracy or all-zero generators do not have any effect here
    
    % take random values for factors
    factors = -1 + 2*rand(size(G,2),N);
    % sample points
    p = c + G * factors;
    
% sampling of extreme random points
elseif strcmp(type,'extreme')
    
    % consider degenerate case
    if rank(G) < n
        Z = Z + (-c);
        p = c * zeros(1,N);
        [S,V,~] = svd([-G,G]);
        d = diag(V);
        ind = find(d > eps);
        if isempty(ind)
            return;
        end
        Z = project(S'*Z,ind);
        temp = randPoint(Z,N,type);
        p(ind,:) = temp;
        p = c + S*p;
        return;
    end
    
    % remove redundant generators
    Z = deleteZeros(Z);
    Z = deleteAligned(Z);
    
    % compute number of zonotope vertices
    q = numberZonoVertices(Z);
    
    if ischar(N) && strcmp(N,'all')
        % return all extreme points
        p = vertices(Z);
        
    elseif 10*N < q
        % generate random vertices
        p = getRandomVertices(Z,N);
        
    elseif N < q
        % select random vertices
        V = vertices(Z);
        ind = randperm(q);
        V = V(:,ind);
        p = V(:,1:N);
        
    else
        % compute vertices and additional points on the boundary
        V = vertices(Z);
        N_ = N - size(V,2);
        V_ = getRandomBoundaryPoints(Z,N_);
        p = [V,V_];
    end
end


end


% Auxiliary Functions -----------------------------------------------------

function V = getRandomVertices(Z,N)
% generate random vertices

    n = dim(Z); m = size(Z.Z,2)-1;
    V = zeros(m,N); cnt = 1; G = generators(Z);
 
    % loop until the desired number of vertices is achieved
    while cnt <= N
        
        % generate random zonotope face
        temp = randperm(m);
        ind = temp(1:n-1);
        Q = G(:,ind);
        c = ndimCross(Q);
        v = sign(c'*G)';
        
        % generate random vertex on the zonotope face
        while true
           v_ = v;
           v_(ind) = sign((-1 + 2*rand(n-1,1)));
           if ~ismember(v_',V','rows')
              V(:,cnt) = v_;
              cnt = cnt + 1;
              break;
           end
        end
    end
    
    % compute vertices
    V = center(Z) + G*V;
end

function V = getRandomBoundaryPoints(Z,N)
% generate random points on the zonotope vertices

    n = dim(Z); m = size(Z.Z,2)-1;
    V = zeros(m,N); G = generators(Z);
 
    % loop until the desired number of vertices is achieved
    for i = 1:N
        
        % generate random zonotope face
        temp = randperm(m);
        ind = temp(1:n-1);
        Q = G(:,ind);
        c = ndimCross(Q);
        r = rand();
        if r > 0.5
           c = -c; 
        end
        V(:,i) = sign(c'*G);
        
        % generate random point on the zonotope face
        V(ind,i) = -1 + 2*rand(n-1,1);
    end
    
    % compute vertices
    V = center(Z) + G*V;
end

function q = numberZonoVertices(Z)
% compute the number of zonotope vertices

    n = dim(Z); m = size(Z.Z,2)-1;
    D = zeros(n,m);
    D(1,:) = 2*ones(1,size(D,2));
    D(:,1) = 2*ones(size(D,1),1);

    for j = 2:size(D,1)
        for k = 2:size(D,2)
            D(j,k) = D(j,k-1) + D(j-1,k-1);
        end
    end
    
    q = D(end,end); 
end

%------------- END OF CODE --------------