function pZ = reduceAdaptive(pZ,diagpercent)
% reduceAdaptive - reduces the zonotope order until a maximum amount of
%    over-approximation defined by the Hausdorff distance between the
%    original zonotope and the reduced zonotope; based on [Thm 3.2,1]
%
% Syntax:
%    pZ = reduceAdaptive(pZ,diagpercent)
%
% Inputs:
%    pZ - polyZonotope object
%    diagpercent - percentage of diagonal of box over-approximation of
%                  polyZonotope (used to compute dHmax)
%
% Outputs:
%    pZ - reduced polyZonotope object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1 0.02 0.003; 0 2 1 0.01 -0.001];
%    GI = [0;0.5];
%    E = [1 0 3 0 1;0 1 1 2 1];
%    pZ = polyZonotope(c,G,GI,E);
%    pZ = reduce(pZ,'adaptive',0.05);
%
% References:
%    [1] Wetzlinger et al. "Adaptive Parameter Tuning for Reachability 
%        Analysis of Nonlinear Systems", HSCC 2021             
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: reduce

% Authors:       Mark Wetzlinger
% Written:       01-October-2020
% Last update:   16-June-2021 (restructure input/output args)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read data from pZ
G = pZ.G;
GI = pZ.GI;
E = pZ.E;
halfs = ~any(mod(E,2),1); % dep. gens with only-even exponents
if any(halfs)
    G(:,halfs) = G(:,halfs) * 0.5; % decrease G if only-even exponents
end
[n,nrG] = size(G);
nrGI = size(GI,2);
Gabs = abs(G);
GIabs = abs(GI);

% set dHmax percentage of diagonal of box(Z)
Gbox = sum([Gabs,GIabs],2); % same as: rad(interval(pZ))
dHmax = (diagpercent * 0.5) * sqrt(sum(Gbox.^2));

m_dep = []; hadapG = [];
if ~isempty(Gabs)
    % order generators for selection by extended girard metrics
    % dep. gens: similar to vanilla reduce
    m_dep = sum(Gabs,1);
    [~,idxDep] = mink(m_dep,nrG);
    % compute dH over-approximation for dep. gens
    % use pZ.G to account for center shift
    % (we use max distance between origin and symm.int.hull)
    % multiplicative factor 2 needed due to deletion of dep. gen from G
    hadapG = 0.5 * vecnorm(cumsum(abs(pZ.G(:,idxDep)),2),2);
    % disect individual cumulative sum to unify below to combined metrics
    hadapG = [hadapG(1),diff(hadapG)];
end

m_indep = []; hadapGI = [];
if ~isempty(GIabs)
    % indep. gens: same as for zonotopes
    [norminf,maxidx] = max(GIabs,[],1); % maxidx for mu vector
    normsum = sum(GIabs,1);
    m_indep = normsum - norminf;
    [~,idxIndep] = mink(m_indep,nrGI);
    % compute additional vector for indep. gens, use linear indexing
    muGIabs = zeros(n,nrGI);
    muGIabs(n*(0:nrGI-1)+maxidx) = norminf;
    % compute new over-approximation of dH
    GIdiag = cumsum(GIabs(:,idxIndep)-muGIabs(:,idxIndep),2);
    hadapGI = 2 * vecnorm(GIdiag,2);
    % disect individual cumulative sus to unify below to combined metrics
    hadapGI = [hadapGI(1),diff(hadapGI)];
end

% order for both together
[~,idxall] = mink([m_dep,m_indep],nrG + nrGI);

hext = zeros(1,nrG + nrGI);
hext(idxall > nrG) = hadapGI;
hext(idxall <= nrG) = hadapG;
hext = cumsum(hext);

% reduce until over-approximation of dH hits dHmax
redUntil = find(hext <= dHmax,1,'last');
idxDep = idxall(idxall(1:redUntil) <= nrG);
idxIndep = idxall(idxall(1:redUntil) > nrG) - nrG;

% delete converted generators
pZ.G(:,idxDep) = [];
E(:,idxDep) = [];
GI(:,idxIndep) = [];
if isempty(idxDep) && isempty(idxIndep)
    pZ.GI = GI;
else
    pZ.GI = [GI, diag(sum([Gabs(:,idxDep),GIabs(:,idxIndep)],2))];
end

% shift pZ.c by center of zonotope converted from dep. gens
if any(halfs) && ~isempty(idxDep)
    temp = find(halfs);
%     if isscalar(temp) 
%         temp = idxDep(idxDep == temp);
%     else
        temp = temp(ismember(temp,idxDep));
%     end
%     if ~isempty(temp)
        pZ.c = pZ.c + sum(G(:,temp),2); % 0.5 factor already done above
%     end
end

% remove all unused dependent factors (empty rows in E)
temp = any(E,2);
pZ.E = E(temp,:);
pZ.id = pZ.id(temp);

end

% ------------------------------ END OF CODE ------------------------------
