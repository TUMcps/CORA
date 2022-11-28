function p = randPoint(pZ,varargin)
% randPoint - generates a random point within a polynomial zonotope
%
% Syntax:  
%    p = randPoint(pZ)
%    p = randPoint(pZ,N)
%    p = randPoint(pZ,N,type)
%    p = randPoint(pZ,'all','extreme')
%
% Inputs:
%    pZ - polyZonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    pZ = polyZonotope([0;0], [2 0 1;1 2 1],[],[1 0 1;0 1 3]);
%    p = randPoint(pZ);
%
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_randPoint('polyZonotope',pZ,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of vars
    p = vars{1}; return
else
    % assign variables
    pZ = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
    % sampling with 'gaussian' is done in contSet method
    if strcmp(type,'gaussian')
        p = randPoint@contSet(pZ,N,type,pr); return
    end
end

% get object properties
m = length(pZ.id); q = size(pZ.Grest,2); 

% compute random points for factor domain interval \alpha \in [-1,1]
dom = interval(-ones(m+q,1),ones(m+q,1));

fac = randPoint(dom,N,type);

% center
p = pZ.c * ones(1,size(fac,2));

% Part 1: dependent generators
if ~isempty(pZ.G)
    for i = 1:size(fac,2)
        p(:,i) = p(:,i) + pZ.G * prod(fac(1:m,i).^pZ.expMat,1)';
    end
end

% Part 2: independent generators
if ~isempty(pZ.Grest)
    p = p + pZ.Grest * fac(m+1:end,:);
end

%------------- END OF CODE --------------