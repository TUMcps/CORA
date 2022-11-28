function p = randPoint(P,varargin)
% randPoint - generates a random point within a polytope
%
% Syntax:  
%    p = randPoint(P)
%    p = randPoint(P,N)
%    p = randPoint(P,N,type)
%    p = randPoint(P,'all','extreme')
%    p = randPoint(P,N,'gaussian',pr)
%
% Inputs:
%    P - mptPolytope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme', or 'gaussian')
%    pr - probability that a value is within the set (only type = 'gaussian')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    P = mptPolytope.generateRandom('Dimension',2);
%    points = randPoint(P,100);
%
%    figure; hold on;
%    plot(P,[1,2],'r');
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
%               19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:---

%------------- BEGIN CODE --------------

    % pre-processing
    [res,vars] = pre_randPoint('mptPolytope',P,varargin{:});
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of vars
        p = vars{1}; return
    else
        % assign variables
        P = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
        % sampling with 'gaussian' is done in contSet method
        if strcmp(type,'gaussian')
            p = randPoint@contSet(P,N,type,pr); return
        end
    end
    
    
    % return all extreme points 
    if ischar(N) && strcmp(N,'all')
        p = vertices(P); return;
    end
    
    % generate random points
    p = zeros(dim(P),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            p(:,i) = randPointNormal(P);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = randPointExtreme(P);
        end
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