function p = randPoint_(C,N,type,varargin)
% randPoint_ - Returns a random point of a capsule
%
% Syntax:
%    p = randPoint_(C)
%    p = randPoint_(C,N)
%    p = randPoint_(C,N,type)
%
% Inputs:
%    C - capsule object
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
% See also: contSet/randPoint, zonotope/randPoint_

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   19-August-2022 (MW, integrate preprocessing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 'all' vertices not supported
if ischar(N) && strcmp(N,'all')
    throw(CORAerror('CORA:notSupported',...
        "Number of vertices 'all' is not supported for class capsule."));
end

% dimension of capsule
n = dim(C);

% initialize points
p = zeros(n,N);

% generate different types of points
if strcmp(type,'standard')
    
    for i = 1:N
        % sample a direction
        dir = -1 + 2*rand(n,1);
        % normalize the direction
        dir = dir / norm(dir,2);
        % randomly choose values for the generator and the radius
        p(:,i) = C.c + (-1 + 2*rand)*C.g + dir * (rand*C.r); 
    end
    
elseif strcmp(type,'extreme')
    
    for i = 1:N
        % sample a direction
        dir = -1 + 2*rand(n,1);
        % normalize the direction
        dir = dir / norm(dir,2);
        % evaluate the support function in that direction
        [~,x] = supportFunc_(C,dir,'upper');
        % choose the point where the support function attains its maximum
        p(:,i) = x; 
    end
    
else
    throw(CORAerror('CORA:noSpecificAlg',type,C));
end

% ------------------------------ END OF CODE ------------------------------
