function p = randPoint(obj,varargin)
% randPoint - Returns a random point of a capsule
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%
% Inputs:
%    obj - capsule object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point inside the capsule
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    p = randPoint(C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
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

    % check for empty set
    if isempty(obj)
        [msg,id] = errEmptySet();
        error(id,msg);
    end
    
    % get object properties
    n = dim(obj);

    % generate different types of extreme points
    if strcmp(type,'standard')
        
        p = zeros(n,N);
        
        for i = 1:N
        
            dir = -1 + 2*rand(n,1);
            dir = dir / norm(dir,2);

            p(:,i) = obj.c + (-1 + 2*rand())*obj.g + dir * (rand()*obj.r); 
        end
        
    elseif strcmp(type,'extreme')
        
        p = zeros(n,N);
        
        for i = 1:N
        
            dir = -1 + 2*rand(n,1);
            dir = dir / norm(dir,2);

            [~,x] = supportFunc(obj,dir);
            
            p(:,i) = x; 
        end
        
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end
end

%------------- END OF CODE --------------