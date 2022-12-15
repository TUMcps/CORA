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

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      23-March-2018
% Last update:  19-August-2022 (MW, integrate standardized pre-processing)
%               09-December-2022 (TL: vectorized computation)
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
n = dim(pZ);
m = length(pZ.id); 
q = size(pZ.Grest,2); 

% compute random points for factor domain interval \alpha \in [-1,1]
dom = interval(-ones(m+q,1),ones(m+q,1));

fac = randPoint(dom,N,type);
N = size(fac, 2); % in case N = 'all'

% center
p = pZ.c;

% dependent generators
if ~isempty(pZ.G)
    fac_dep = fac(1:m, :);
    fac_dep = reshape(fac_dep, m, 1, N);
    fac_dep = prod(fac_dep .^ pZ.expMat, 1);
    dep = sum(pZ.G .* fac_dep, 2);
    dep = reshape(dep, n, N);
    p = p + dep;
end

% independent generators
if ~isempty(pZ.Grest)
    fac_ind = fac(m+1:end, :);
    fac_ind = reshape(fac_ind, 1, [], N);
    ind = sum(pZ.Grest .* fac_ind, 2);
    ind = reshape(ind, n, N);
    p = p + ind;
end

%------------- END OF CODE --------------