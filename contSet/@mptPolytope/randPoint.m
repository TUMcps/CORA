function p = randPoint(obj,varargin)
% randPoint - generates a random point within a polytope
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
%    obj - mptPolytope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', or 'gaussian')
%    pr - probability that a value is within the set (only type = 'gaussian')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    poly = mptPolytope.generateRandom(2);
%
%    points = randPoint(poly,100);
%
%    figure; hold on;
%    plot(poly,[1,2],'r');
%    plot(points(1,:),points(2,:),'.k','MarkerSize',10);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  25-June-2021 (MP, add type gaussian)
% Last revision:---

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
                    if nargin > 3  && strcmp(type,'gaussian')
                        pr = varargin{3};
                    end
                end
            end
            
        elseif ischar(varargin{1}) && any(strcmp(varargin{1},types))
            % second input argument is type (sampling method)
            type = varargin{1};
            if nargin > 2 && strcmp(type,'gaussian')
                pr = varargin{2};
            end
            
        else
            [msg,id] = errWrongInput('N or type');
            error(id,msg);
        end
    end
    
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(obj); return;
    end
    
    % generate random points
    p = zeros(dim(obj),N);
    
    
    if strcmp(type,'gaussian')
        if nargin == 3
            p = randPoint@contSet(obj,type,pr);
        else
            p = randPoint@contSet(obj,N,type,pr);
        end
    elseif strcmp(type,'standard')
        for i = 1:N
            p(:,i) = randPointNormal(obj);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = randPointExtreme(obj);
        end
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end
end


% Auxiliary Functions -----------------------------------------------------

function p = randPointNormal(P)
% generate random point within the polytope

    % draw n+1 random extreme points
    n = dim(P);
    points = zeros(n,n+1);
    
    for i = 1:size(points,2)
       points(:,i) = randPointExtreme(P); 
    end

    % interpolate betwenn the points
    delta = rand(n+1,1);
    delta = delta./sum(delta);
    
    p = points*delta;
end

function p = randPointExtreme(P)
% generate random point on the polytope boundary
        
    % check if vertex representation is available
    if P.P.hasVRep
        
       % select random vertex
       V = P.P.V;
       ind = randi([1,size(V,1)]);
       p = V(ind,:)';
       
    else
       
       % center polytope at origin
       c = center(P);
       P = P - c;
       
       % select random direction
       n = length(c);
       d = rand(n,1) - 0.5*ones(n,1);
       d = d./norm(d);
        
       % compute farest point in this direction that is still in polytope
       [~,x] = supportFunc(P,d);
       p = x + c;
    end
end

%------------- END OF CODE --------------