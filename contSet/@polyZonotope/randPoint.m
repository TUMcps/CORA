function p = randPoint(obj,varargin)
% randPoint - generates a random point within a polynomial zonotope
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%
% Inputs:
%    obj - polyZonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    pZ = polyZonotope([0;0], [2 0 1;1 2 1],[],[1 0 1;0 1 3]);
%    
%    p = randPoint(pZ);
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    if nargin > 1 && ~isempty(varargin{1})
       N = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
       type = varargin{2}; 
    end

    % get object properties
    m = length(obj.id); q = size(obj.Grest,2); 
    
    % compute random points for factor domain interval \alpha \in [-1,1]
    dom = interval(-ones(m+q,1),ones(m+q,1));

    fac = randPoint(dom,N,type);

    % center
    p = obj.c * ones(1,size(fac,2));

    % Part 1: dependent generators
    if ~isempty(obj.G)
        for i = 1:size(fac,2)
            p(:,i) = p(:,i) + obj.G * prod(fac(1:m,i).^obj.expMat,1)';
        end
    end

    % Part 2: independent generators
    if ~isempty(obj.Grest)
        p = p + obj.Grest * fac(m+1:end,:);
    end
end

%------------- END OF CODE --------------