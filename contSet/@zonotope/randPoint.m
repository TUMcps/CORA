function p = randPoint(obj,varargin)
% randPoint - generates random points within a zonotope
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,type)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%    p = randPoint(obj,N,'gaussian',pr)
%    p = randPoint(obj,'gaussian',pr)
%
% Inputs:
%    obj - zonotope object
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
% Last revision: ---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    types = {'standard','extreme','gaussian'};
    
    if nargin > 1 && ~isempty(varargin{1})
        if (isnumeric(varargin{1}) && isscalar(varargin{1})) || ...
                (ischar(varargin{1}) && strcmp(varargin{1},'all'))
            % second input argument is number of points
            N = varargin{1};
            
            if nargin > 2
                if ischar(varargin{2}) && any(strcmp(varargin{2},types))
                    type = varargin{2};
                    if strcmp(type,'gaussian')
                    	pr = varargin{3};
                    end
                end
            end
            
        elseif ischar(varargin{1}) && any(strcmp(varargin{1},types))
            % second input argument is type (sampling method)
            type = varargin{1};
            if strcmp(type,'gaussian')
                pr = varargin{2};
            end
            
        else
            [msg,id] = errWrongInput('N or type');
            error(id,msg);
        end
    end
   
    
    % generate different types of extreme points
    if strcmp(type,'gaussian')
        if nargin == 3
            p = randPoint@contSet(obj,type,pr);
        else
            p = randPoint@contSet(obj,N,type,pr);
        end
    else
        % get object properties
        c = center(obj); G = generators(obj); n = dim(obj);
        
        % empty set
        if n == 0
            [msg,id] = errEmptySet();
            error(id,msg);
        end
        
        if strcmp(type,'standard')
            
            factors = -1 + 2*rand(size(G,2),N);
            p = c + G * factors;
            
        elseif strcmp(type,'extreme')
            
            % consider degenerate case
            if rank(G) < n
                obj = obj + (-c);
                p = c * zeros(1,N);
                [S,V,~] = svd([-G,G]);
                d = diag(V);
                ind = find(d > eps);
                if isempty(ind)
                    return;
                end
                obj = project(S'*obj,ind);
                temp = randPoint(obj,N,type);
                p(ind,:) = temp;
                p = c + S*p;
                return;
            end
            
            % remove redundant generators
            obj = deleteZeros(obj);
            obj = deleteAligned(obj);
            
            % compute number of zonotope vertices
            q = numberZonoVertices(obj);
            
            % return all extreme point
            if ischar(N) && strcmp(N,'all')
                
                p = vertices(obj);
                
                % generate random vertices
            elseif 10*N < q
                
                p = getRandomVertices(obj,N);
                
                % select random vertices
            elseif N < q
                
                V = vertices(obj);
                ind = randperm(q);
                V = V(:,ind);
                p = V(:,1:N);
                
                % compute vertices and additional points on the boundary
            else
                
                V = vertices(obj);
                N_ = N - size(V,2);
                V_ = getRandomBoundaryPoints(obj,N_);
                p = [V,V_];
            end
        else
            [msg,id] = errWrongInput('type');
            error(id,msg);
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