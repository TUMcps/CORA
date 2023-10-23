function [Z,dHerror,gredIdx] = reduceAdaptive(Z,diagpercent,varargin)
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
%               zonotope (used to compute dHmax) [0,1]
%
% Outputs:
%    Z - reduced zonotope
%    dHerror - Hausdorff distance between Z and reduced Z
%    gredIdx - index of reduced generators
%
% Example: 
%    Z = zonotope([1;0],[1 3 2 -1 0.03 0.02 -0.1; 2 0 -1 1 0.02 -0.01 0.2]);
%    Z = reduce(Z,'adaptive',10);
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

% default type
type = 'girard';
if nargin == 3
    if ischar(varargin{1}) && any(ismember(varargin{1},{'penven','girard'}))
        type = varargin{1};
    end
end

dHerror = 0;
G = generators(Z);
if isempty(G)
    dHerror = 0;
    gredIdx = [];
    return;
end
Gabs = abs(G);

% compute maximum admissible dH
Gbox = sum(Gabs,2);
dHmax = (diagpercent * 2) * sqrt(sum(Gbox.^2));

[n,nrG] = size(G);

if strcmp(type, 'penven')
    % dummy value (plz fix Adri)
    gredIdx = 42;
    
    % Compute spherical (naive) bound for each generator
    naive = sqrt(sum(Gabs.^2,1));
    % Compute Le Penven bound for each generator, by multiplying the last bit to the spherical bound
    penven = naive .* sqrt(2*(abs(1-  sum(Gabs.^4,1)./sum(Gabs.^2,1).^2 )  ));
    % Compare them
    resulting = min([naive;penven]);

    % Sort them
    [h,idx] = mink(resulting,nrG);

    if ~any(h)
        % no generators or all are h=0
        newG = diag(Gbox);
        % Z.c = Z.c;
        Z.G = newG(:,any(newG,1));
        return
    end

    % box generators with h = 0
    hzeroIdx = idx(h==0);
    Gzeros = sum(Gabs(:,hzeroIdx),2);
    last0Idx = numel(hzeroIdx);
    gensred = Gabs(:,idx(last0Idx+1:end));

    % Take cumsum of the bounds for each generator, take the first few
    s = cumsum(h);
    redIdx = find(s(last0Idx+1:end) <= dHmax, 1, 'last');
    if isempty(redIdx)
        redIdx = 0;
        dHerror = 0;
        gredIdx = hzeroIdx;
    else
        dHerror = h(redIdx);
        gredIdx = idx(1:length(hzeroIdx)+redIdx);
    end
    
else % 'girard'
    % select generators using 'girard'
    norminf = max(Gabs,[],1);               % faster than: vecnorm(G,Inf);
    normsum = sum(Gabs,1);                  % faster than: vecnorm(G,1);
    [h,idx] = mink(normsum - norminf,nrG);

    if ~any(h)
        % no generators or all are h=0
        newG = diag(Gbox);
        % Z.c = Z.c
        Z.G = newG(:,any(newG,1));
        % no over-approximation
        dHerror = 0;
        gredIdx = idx;
        return
    end

    % box generators with h = 0
    hzeroIdx = idx(h==0);
    Gzeros = sum(Gabs(:,hzeroIdx),2);
    last0Idx = numel(hzeroIdx);
    gensred = Gabs(:,idx(last0Idx+1:end));

    [maxval,maxidx] = max(gensred,[],1);
    % use linear indexing
    mugensred = zeros(n,nrG-last0Idx);
    cols = n*(0:nrG-last0Idx-1);
    mugensred(cols+maxidx) = maxval;
    % compute new over-approximation of dH
    gensdiag = cumsum(gensred-mugensred,2);
    h = 2 * vecnorm(gensdiag,2); %sqrt(sum(gensdiag.^2,1))
    % index until which gens are reduced
    redIdx = find(h <= dHmax,1,'last');
    if isempty(redIdx)
        redIdx = 0;
        dHerror = 0;
        gredIdx = hzeroIdx;
    else
        dHerror = h(redIdx);
        gredIdx = idx(1:length(hzeroIdx)+redIdx);
    end
end


Gred = sum(gensred(:,1:redIdx),2);
% prior version: no sorting
% Gunred = G(:,idx(last0Idx+redIdx+1:end));
% sort to keep correspondances!
Gunred = G(:,sort(idx(last0Idx+redIdx+1:end)));

% Z.c = Z.c;
Z.G = [Gunred,diag(Gred+Gzeros)];

% just for performance evaluation -----------------------------------------

% compute actual dH
% gred = G(:,idx(1:last0Idx+redIdx));
% Ztest = zonotope([zeros(n,1),gred]);
% shrinkFactors = 1 ./ sum(generators(box(Ztest)),2);
% Ztest = enlarge(Ztest,shrinkFactors);
% [realdH,hcomp] = hausdorffBox(Ztest);

end

% ------------------------------ END OF CODE ------------------------------
