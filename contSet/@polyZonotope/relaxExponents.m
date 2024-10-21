function [pZ_relax,dig] = relaxExponents(pZ,varargin)
% relaxExponents - relaxes the exponents of a polynomial zonotope
%
% Syntax:
%    pZ = relaxExponents(pZ,eta,method)
%
% Inputs:
%    pZ - polyZonotope object
%    eta - exponent threshold
%    method - 'greedy' or 'all'
%
% Outputs:
%    pZ_relax - relaxed polyZonotope object
%    dig - relaxation graph
%
% References:
%    [1] Ladner, T., et al. (2024). Exponent relaxation of polynomial zonotopes 
%        and its applications in formal neural network verification. AAAI.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       03-July-2023
% Last update:   ---
% Last revision: 03-July-2023 (combined functions)

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,3);

[eta,method] = setDefaultValues({1,'greedy'},varargin);
inputArgsCheck({ ...
    {pZ,'att','polyZonotope'}, ...
    {eta,'att','numeric',{'scalar','integer','nonnegative'}}, ...
    {method,'str',{'greedy','all'}}, ...
})

% construct graph
dig = aux_constructGraph(pZ,eta,method);

if strcmp(method,'all')
    % only return graph
    pZ_relax = [];
    return
end

% merge generators according to graph
pZ_relax = aux_mergeGenerators(pZ,dig);

if nargout < 2
    clear dig;
end

end


% Auxiliary functions -----------------------------------------------------

function dig = aux_constructGraph(pZ,eta,method)

% read exponent matrix
n = dim(pZ);
G = pZ.G;
E = pZ.E;
[~,h] = size(E);

% store error bounds
maxE = max(E,[],"all");
errorl = zeros(maxE+1,maxE+1);
erroru = zeros(maxE+1,maxE+1);

% successor vector and edge weights (type='greedy')
S = zeros(1,h);
W = zeros(1,h);

% adjacency matrix (type='all')
A = sparse(h,h);

H = 1:h;
for i=H
    % compute generators with at most one difference to generator i
    idxMax1Diff = sum((E(:,i) ~= E),1) == 1;

    % compute similar ones
    diff = E(:,idxMax1Diff) - E(:,i);
    idx_sim = all(diff >= 0,1) ... % from high to low
        & all((E(:,idxMax1Diff) >= eta & E(:,i) >= eta) | diff == 0,1) ... % >= minE
        & (...
            mod(sum(diff,1),2) == 0 ... % both odd or both even
            | all((abs(diff) == 1 & E(:,i) == 0) | diff == 0,1) ... % a^1 -> a^0
        );

    % update indices
    idxMax1Diff(idxMax1Diff) = idx_sim;

    % add edges (j -> i) to graph
    for j=H(idxMax1Diff)
        
        if n == 1
            gi = G(:,i);
            gj = G(:,j);
            if sign(gi) == sign(gj)
                continue
            end
        end

        % read exponents
        idx = E(:,i) ~= E(:,j);
        ei = E(idx,i);
        ej = E(idx,j);
            
        % compute approx error
        dl = errorl(ej+1,ei+1); % reuse error if present
        du = erroru(ej+1,ei+1); % reuse error if present
        if du-dl == 0
            % recompute error
            pi = [unitvector(1,ei)',0];
            pj = [unitvector(1,ej)',0];
            [dl,du] = minMaxDiffPoly(pi,pj,-1,1);
            % save error for reusage
            errorl(ej+1,ei+1) = dl;
            erroru(ej+1,ei+1) = du;
        end
        
        dc = (dl+du)/2;
        if all(mod(E(:,i),2) == 0)
            % all even exponents, error can be halfed
            % shift center
            % c = c + dc*G(:,j);
            % error is radius then
            di = du-dc; 
    
        elseif abs(dc) < 1e-8 && (ei ~= 0)
            % error is symmetric
            di = (du-dl)/2;
        else
            % over-approximate
            di = du-dl;
        end

        % consider direction using scalar product
        % di = di + (1+(gi'*gj)/(vecnorm(gi)*vecnorm(gj)));

        if strcmp(method,'greedy')
            % find current successor of generator j (j -> k)
            k = S(j);
            dk = W(j);

            if k == 0 % no successor
                % add edge (j -> i)
                S(j) = i;
                W(j) = di;
    
            else
                % compare successor's edge weight with current
                if dk > di                    
                    % change (j -> k) to (j -> i)
                    S(j) = i;
                    W(j) = di;
                end
            end
        elseif strcmp(method,'all')
            A(i,j) = di;
            A(j,i) = di;
        end
    end
end

% init graph
if strcmp(method,'greedy')
    % init directed graph
    idx = S>0;
    dig = digraph(H(idx),S(idx),W(idx));
elseif strcmp(method,'all')
    % init undirected graph
    dig = graph(A);
end

end

% ---

function pZ_relax = aux_mergeGenerators(pZ,dig)

if height(dig.Edges) == 0
    % nothing to merge
    pZ_relax = pZ;
    return;
end

% read properties
c = pZ.c;
G = pZ.G;
GI = pZ.GI;
E = pZ.E;
id = pZ.id;

% store error bounds
maxE = max(E,[],"all");
errorl = zeros(maxE+1,maxE+1);
erroru = zeros(maxE+1,maxE+1);

% preallocate GI
n = length(c);
[n_I,q] = size(GI);
GI_ext = zeros(n,q+height(dig.Edges));
GI_ext(1:n_I,1:q) = GI;

% determine generator merging order 
mergeOrder = zeros(1,height(dig.Nodes));
cnt = 0;

% using inverse bfs
dig = flipedge(dig);
% find root generators, which are now leafs
leafs = find(indegree(dig) == 0 & outdegree(dig) > 0)';
for i = leafs
    % breadth-first search
    nodes = bfsearch(dig,i)';

    % reverse order
    nodes = fliplr(nodes); 

    % remove root gnerator
    nodes(end) = [];

    % add to nodeOrder
    mergeOrder(cnt+1:cnt+length(nodes)) = nodes;
    cnt = cnt + length(nodes);
end
dig = flipedge(dig);
mergeOrder(cnt+1:end) = [];

% process nodes
for i=mergeOrder
    % read generators and exponents
    Gi = G(:,i);
    Ei = E(:,i);

    % find successor
    j = successors(dig,i);
    Ej = E(:,j);

    % compute approx error
    idx = Ei ~= Ej;        
    if nnz(idx) ~= 1
        throw(CORAerror("CORA:wrongValue",'second','Exponents of generators should only differ at one dependent factor if edge in graph is present.'))
    end

    % read exponents
    ei = Ei(idx);
    ej = Ej(idx);
    
    % compute approx error
    dl = errorl(ej+1,ei+1); % reuse error if present
    du = erroru(ej+1,ei+1); % reuse error if present
    if du-dl == 0
        % recompute error
        pi = [unitvector(1,ei)',0];
        pj = [unitvector(1,ej)',0];
        [dl,du] = minMaxDiffPoly(pi,pj,-1,1);
        % save error for reusage
        errorl(ej+1,ei+1) = dl;
        erroru(ej+1,ei+1) = du;
    end
    
    dc = (dl+du)/2;
    if all(mod(Ej,2) == 0)
        % all even exponents, error can be halfed
        % shift center
        c = c + dc*Gi;
        % error is radius then
        d = du-dc; 

    elseif abs(dc) < 1e-8 && (ej ~= 0)
        % error is symmetric
        d = (du-dl)/2;
    else
        % over-approximate
        d = du-dl;
    end

    % read generator j
    Gj = G(:,j);
    
    % add generator i to generator j 
    G(:,j) = Gj + Gi;
    % add approx error
    q = q+1;
    GI_ext(:,q) = d*Gi;
    
    % 'delete' generator i (but preserve indices)
    G(:,i) = Gi * 0;
    E(:,i) = Ei * 0;
end
    
% remove 'deleted' all-zero generators
idx = all(G == 0,1);
G = G(:,~idx);
E = E(:,~idx);

% construct resulting polyZonotope
pZ_relax = polyZonotope(c,G,GI_ext,E,id);

end

% ------------------------------ END OF CODE ------------------------------
