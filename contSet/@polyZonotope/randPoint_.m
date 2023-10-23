function p = randPoint_(pZ,N,type,varargin)
% randPoint_ - generates a random point within a polynomial zonotope
%
% Syntax:
%    p = randPoint_(pZ)
%    p = randPoint_(pZ,N)
%    p = randPoint_(pZ,N,type)
%    p = randPoint_(pZ,'all','extreme')
%
% Inputs:
%    pZ - polyZonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%           Note that for 'extreme', the generated points might not all be
%           extremal.
%
% Outputs:
%    p - random points in R^n
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
% See also: contSet/randPoint, zonotope/randPoint_

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       23-March-2018
% Last update:   19-August-2022 (MW, integrate standardized pre-processing)
%                09-December-2022 (TL, vectorized computation)
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

if ~(strcmp(type,'standard') || strcmp(type,'extreme'))
    throw(CORAerror('CORA:noSpecificAlg',type,pZ));
end

% get object properties
n = dim(pZ);
m = length(pZ.id); 
q = size(pZ.GI,2); 

% compute random points for factor domain interval \alpha \in [-1,1]
dom = interval(-ones(m+q,1),ones(m+q,1));

fac = randPoint_(dom,N,type);
N = size(fac, 2); % in case N = 'all'

% center
p = pZ.c;

% dependent generators
if ~isempty(pZ.G)
    fac_dep = fac(1:m, :);
    fac_dep = reshape(fac_dep, m, 1, N);    
    fac_dep = prod(fac_dep .^ pZ.E, 1);
    dep = sum(pZ.G .* fac_dep, 2);
    dep = reshape(dep, n, N);
    p = p + dep;
end

% independent generators
if ~isempty(pZ.GI)
    fac_ind = fac(m+1:end, :);
    fac_ind = reshape(fac_ind, 1, [], N);
    ind = sum(pZ.GI .* fac_ind, 2);
    ind = reshape(ind, n, N);
    p = p + ind;
end

% ------------------------------ END OF CODE ------------------------------
