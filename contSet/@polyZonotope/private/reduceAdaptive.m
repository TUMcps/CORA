function pZ = reduceAdaptive(pZ,diagpercent)
% reduceAdaptive - reduces the zonotope order until a maximum amount of
%    over-approximation defined by the Hausdorff distance between the
%    original zonotope and the reduced zonotope; based on [Thm 3.2,1]
%
% Syntax:  
%    Z = reduceAdaptive(Z,diagpercent)
%
% Inputs:
%    Z - zonotope object
%    diagpercent - percentage of diagonal of box over-approximation of
%                  polyZonotope (used to compute dHmax)
%
% Outputs:
%    Z - reduced zonotope
%
% Example: 
%    pZ = polyZonotope.generateRandom(15);
%    pZ = reduce(Z,'adaptive',0.05);
%    
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

% Author:       Mark Wetzlinger
% Written:      01-October-2020
% Last update:  16-June-2021 (restructure input/output args)
% Last revision: ---

%------------- BEGIN CODE --------------

% read data from pZ
G = pZ.G;
Grest = pZ.Grest;
expMat = pZ.expMat;
halfs = ~any(mod(expMat,2),1); % dep. gens with only-even exponents
if any(halfs)
    G(:,halfs) = G(:,halfs) * 0.5; % decrease G if only-even exponents
end
[n,nrG] = size(G);
nrGrest = size(Grest,2);
Gabs = abs(G);
Grestabs = abs(Grest);

% set dHmax percentage of diagonal of box(Z)
Gbox = sum([Gabs,Grestabs],2); % same as: rad(interval(pZ))
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

m_indep = []; hadapGrest = [];
if ~isempty(Grestabs)
    % indep. gens: same as for zonotopes
    [norminf,maxidx] = max(Grestabs,[],1); % maxidx for mu vector
    normsum = sum(Grestabs,1);
    m_indep = normsum - norminf;
    [~,idxIndep] = mink(m_indep,nrGrest);
    % compute additional vector for indep. gens, use linear indexing
    muGrestabs = zeros(n,nrGrest);
    muGrestabs(n*(0:nrGrest-1)+maxidx) = norminf;
    % compute new over-approximation of dH
    Grestdiag = cumsum(Grestabs(:,idxIndep)-muGrestabs(:,idxIndep),2);
    hadapGrest = 2 * vecnorm(Grestdiag,2);
    % disect individual cumulative sus to unify below to combined metrics
    hadapGrest = [hadapGrest(1),diff(hadapGrest)];
end

% order for both together
[~,idxall] = mink([m_dep,m_indep],nrG + nrGrest);

hext = zeros(1,nrG + nrGrest);
hext(idxall > nrG) = hadapGrest;
hext(idxall <= nrG) = hadapG;
hext = cumsum(hext);

% reduce until over-approximation of dH hits dHmax
redUntil = find(hext <= dHmax,1,'last');
idxDep = idxall(idxall(1:redUntil) <= nrG);
idxIndep = idxall(idxall(1:redUntil) > nrG) - nrG;

% delete converted generators
pZ.G(:,idxDep) = [];
expMat(:,idxDep) = [];
Grest(:,idxIndep) = [];
if isempty(idxDep) && isempty(idxIndep)
    pZ.Grest = Grest;
else
    pZ.Grest = [Grest, diag(sum([Gabs(:,idxDep),Grestabs(:,idxIndep)],2))];
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

% remove all unused dependent factors (empty rows in expMat)
temp = any(expMat,2);
pZ.expMat = expMat(temp,:);
pZ.id = pZ.id(temp);

end

