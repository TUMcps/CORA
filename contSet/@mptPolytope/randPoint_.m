function p = randPoint_(P,N,type,varargin)
% randPoint_ - generates a random point within a polytope
%
% Syntax:  
%    p = randPoint_(P)
%    p = randPoint_(P,N)
%    p = randPoint_(P,N,type)
%    p = randPoint_(P,'all','extreme')
%
% Inputs:
%    P - mptPolytope object
%    N - number of random points
%    type - type of the random point ('standard', 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    P = mptPolytope([1 1; -2 1; 0 -1],[1;1;1]);
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
% See also: zonotope/randPoint_

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  25-June-2021 (MP, add type gaussian)
%               19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:27-March-2023 (MW, rename randPoint_)

%------------- BEGIN CODE --------------
    
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
       [~,x] = supportFunc_(P,d,'upper');
       p = x + c;
    end
end

%------------- END OF CODE --------------