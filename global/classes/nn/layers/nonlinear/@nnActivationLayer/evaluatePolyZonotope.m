function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
% evaluatePolyZonotope - evaluates the activation layer on a polyZonotope
%
% Syntax:
%    [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
%
% Inputs:
%    c, G, GI, E, id, id_, ind, ind_ - parameters of polyZonotope
%    evParams - parameter for NN evaluation
%
% Outputs:
%    updated [c, G, GI, E, id, id_, ind, ind_]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnActivationLayer/evaluatePolyZonotopeNeuron

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   05-April-2022
%                23-June-2022 (performance optimizations)
%                05-December-2022 (readability through aux functions)
%                16-February-2023 (combined approx_type)
%                21-March-2023 (bugfix GI)
%                01-August-2023 (simplified order reduction)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init --------------------------------------------------------------------

% init bounds
n = size(G, 1);
if ~all(size(obj.l) == [n, 1])
    obj.l = nan(n, 1);
end
if ~all(size(obj.u) == [n, 1])
    obj.u = nan(n, 1);
end

% prepare properties for refinement
if ~all(size(obj.order) == size(c))
    obj.order = ones(size(c)) .* max(obj.order);
end
obj.do_refinement = ones(size(c));

% due to neuron merging
if n == 0
    return
end

% prepare neuron-wise computation -----------------------------------------

% order reduction prior to the evaluation
[c, G, GI, E, id, id_, ind, ind_] = aux_preOrderReduction(obj, c, G, GI, E, id, id_, ind, ind_, evParams);

% init output sizes using maxOrder
maxOrder = max(obj.order);

if evParams.sort_exponents
    % sort columns of exponent matrix
    [G, E] = aux_sort(obj, G, E, maxOrder);
end

% compute output sizes per order (improves performance)
[G_start, G_end, G_ext_start, G_ext_end] = nnHelper.getOrderIndicesG(G, maxOrder);
[~, GI_end, ~, ~] = nnHelper.getOrderIndicesGI(GI, G, maxOrder);

% preallocate output sizes
c_out = zeros(n, 1);
G_out = zeros(n, G_end(end));
GI_out = zeros(n, GI_end(end));
E_out = aux_computeE_out(obj, E, maxOrder, G_start, G_end);
d = zeros(n, 1);

% loop over all neurons ---------------------------------------------------

for j = 1:n
    evParams.j = j;
    order_j = obj.order(j);

    [c_out(j), G_out_j, GI_out_j, d(j)] = ...
        obj.evaluatePolyZonotopeNeuron( ...
        c(j), G(j, :), GI(j, :), E, E_out, order_j, ...
        ind, ind_, evParams ...
        );

    G_out(j, 1:length(G_out_j)) = G_out_j;
    GI_out(j, 1:length(GI_out_j)) = GI_out_j;
end

if evParams.sort_exponents
    % make sure columns of E_out remain sorted
    [G_out, E_out] = aux_sortPost(obj, G_out, E_out, maxOrder, ...
        G_start, G_end, G_ext_start, G_ext_end);
end

% compute final output ----------------------------------------------------

c = c_out;
G = G_out;
GI = GI_out;
E = E_out;

% order reduction post to the evaluation
[c, G, GI, E, id, d, id_] = aux_postOrderReduction(obj, c, G, GI, E, id, id_, d, evParams);

% add approximation error
[G, GI, E, id, id_] = aux_addApproxError(obj, d, G, GI, E, id, id_, evParams);

% update indices of all-even exponents (for zonotope encl.)
ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
ind_ = setdiff(1:size(E, 2), ind);

end


% Auxiliary functions -----------------------------------------------------

function [c, G, GI, E, id, id_, ind, ind_] = aux_preOrderReduction(obj, c, G, GI, E, id, id_, ind, ind_, evParams)
% order reduction prior to the evaluation

% read number of generators
[n,h] = size(G);
q = size(GI,2);

% read max number of generators
nrMaxGen = evParams.num_generators;
nrMaxGen = min([h+q, nrMaxGen]); % in case evParams.num_generators is empty

if evParams.do_pre_order_reduction && ~isempty(evParams.max_gens_post)
    % consider order reduction
    max_order = max(obj.order);
    nrMaxGenOrderRed = nthroot(evParams.max_gens_post, max_order);
    nrMaxGen = min(nrMaxGen, nrMaxGenOrderRed);
end

if h+q > nrMaxGen
    % reduce
    [c, G, GI, E, id, d] = nnHelper.reducePolyZono(c, G, GI, E, id, nrMaxGen, obj.sensitivity);

    % add to GI
    D = diag(d);
    GI = [GI,D(:,d > 0)];

    % update number of generators
    h = size(G,2);
    q = size(GI,2);

    % update max id
    id_ = max(max(id), id_);
end

% restructure pZ s.t. there remain no independent generators
if size(GI, 2) > 0 && evParams.remove_GI
    % restructure GI to G
    G = [G,GI];
    E = blkdiag(E, eye(q));
    id = [id; id_ + (1:q)'];

    GI = zeros(n, 0);
end

% update auxiliary variables
id_ = max(max(id), id_);
ind = find(prod(ones(size(E))- ...
    mod(E, 2), 1) == 1);
ind_ = setdiff(1:size(E, 2), ind);

end

function [G, E] = aux_sort(obj, G, E, maxOrder)
    % sort columns of exponent matrix in lexigraphic ordering
    if maxOrder > 1
        [E, idx] = sortrows(E');
        E = E';
        G = G(:, idx);
    end
end
 
function [G, E] = aux_sortPost(obj, G, E, maxOrder, G_start, G_end, G_ext_start, G_ext_end)
    % sort columns of exponent matrix in lexigraphic ordering
    if maxOrder > 1
        % sort within PZ^i
        for i=2:maxOrder % skip PZ^1, already sorted
            G_i = G(:, G_start(i):G_end(i));
            E_i = E(:, G_start(i):G_end(i));


            i1 = floor(i/2);
            i2 = ceil(i/2);
            i1len = G_ext_end(i1)-G_ext_start(i1)+1;
            i2len = G_ext_end(i2)-G_ext_start(i2)+1;

            Gs_i = cell(1, i1);
            Es_i = cell(1, i1);
            if i1 == i2
                cnt = 1;
                for j = i1len:-1:1
                    Gs_i{j} = G_i(:, cnt:cnt+j-1);
                    Es_i{j} = E_i(:, cnt:cnt+j-1);
                    cnt = cnt + j;
                end
            else
                idx = num2cell(1:i1len);
                Gs_i = cellfun(@(i) G_i(:, (i-1)*i2len+(1:i2len)), ...
                    idx, 'UniformOutput', false);
                Es_i = cellfun(@(i) E_i(:, (i-1)*i2len+(1:i2len)), ...
                    idx, 'UniformOutput', false);
            end

            [G_i, E_i] = aux_sortPlus(obj, Gs_i, Es_i, size(G_i, 2));
            h_i = size(G_i, 2); % new number of generators
            G(:, G_start(i)-1+(1:h_i)) = G_i;
            G(:, G_start(i)+h_i:G_end(i)) = 0; % set unused to 0
            E(:, G_start(i)-1+(1:h_i)) = E_i;
        end

        % sort between PZ^i
        idx = num2cell(1:maxOrder);
        Gs = cellfun(@(i) G(:, G_start(i):G_end(i)), ...
            idx, 'UniformOutput', false);
        Es = cellfun(@(i) E(:, G_start(i):G_end(i)), ...
            idx, 'UniformOutput', false);
        [G, E] = aux_sortPlus(obj, Gs, Es, size(G, 2));
    end
end
 
function [G, E] = aux_sortPlus(obj, Gs, Es, h)
    % sorts a list of sorted exponent matrices
    if isempty(Gs)
        G = []; E = [];
        return
    end

    % init heap
    heapInit = cell(1, length(Es));
    for i=1:length(Es)
        minObj = struct;
        minObj.idx = 1;
        minObj.i = i;
        minObj.key = Es{minObj.i}(:, minObj.idx)'; % key is row vector
        heapInit{i} = minObj;
    end

    G = zeros(size(Gs{1}, 1), h);
    E = zeros(size(Es{1}, 1), h);
    idx = 1;

    heapObj = nnHelper.heap(heapInit);
    while ~isempty(heapObj)
        minObj = min(heapObj);
        G_i = Gs{minObj.i};
        G(:, idx) = G_i(:, minObj.idx);
        E(:, idx) = minObj.key';
        idx = idx+1;

        if(minObj.idx+1 <= size(G_i, 2))
            % check if last entry of current min Gs/Es is smaller than
            % current min of heap
            E_i = Es{minObj.i};
            minObjNext = min(heapObj);
            [~, ids] = sortrows([minObjNext.key; E_i(:, end)']);
            if ids(1) == 1
                % add next generator
                minObjNew = struct;
                minObjNew.idx = minObj.idx+1;
                minObjNew.i = minObj.i;
                minObjNew.key = Es{minObjNew.i}(:, minObjNew.idx)';
                replaceMin(heapObj, minObjNew);
            else
                % add all remaining generators
                remGens = size(G_i, 2)-minObj.idx;
                G(:, idx:idx+remGens-1) = G_i(:, end-remGens+1:end);
                E(:, idx:idx+remGens-1) = E_i(:, end-remGens+1:end);
                idx = idx + remGens;
                pop(heapObj);
            end
        else
            % remove current min
            pop(heapObj);
        end
    end
            
    [G, E] = aux_merge(obj, G, E);
end

function [G, E] = aux_merge(obj, G, E)
    % merge a sorted, non-regular exponent matrix
    i_out = 1;
    h = size(G, 2);

    for i=2:h
        if E(:, i_out) == E(:, i)
            G(:, i_out) = G(:, i_out) + G(:, i);
        else
            i_out = i_out + 1;
            if i_out ~= i
                G(:, i_out) = G(:, i);
                E(:, i_out) = E(:, i);
            end
        end
    end

    G(:, i_out+1:end) = [];
    E(:, i_out+1:end) = [];
end

function E_out = aux_computeE_out(obj, E, order, G_start, G_end)
% comput output exponential matrix

% init
E_out = zeros(size(E, 1), G_end(end));
E_ext = cell(1, order);
E_out(:, 1:G_end(1)) = E;
E_ext{1} = E;

for i = 2:order
    % Note that e.g., G2 * G3 = G5 -> E2 + E3 = E5
    i1 = floor(i/2);
    i2 = ceil(i/2);

    Ei1_ext = E_ext{i1};
    Ei2_ext = E_ext{i2};
    Ei = nnHelper.calcSquaredE(Ei1_ext, Ei2_ext, i1 == i2);
    Ei_ext = [Ei1_ext, Ei2_ext, Ei];

    if i1 < i2
        % free memory
        E_ext{i1} = [];
    end

    E_out(:, G_start(i):G_end(i)) = Ei;
    E_ext{i} = Ei_ext;
end
end

function [c, G, GI, E, id, d, id_] = aux_postOrderReduction(obj, c, G, GI, E, id, id_, d, evParams)
% order reduction post to the evaluation

% read max number of generators
nrGen = evParams.num_generators;

if ~isempty(nrGen) && nrGen < size(G, 2) + size(GI, 2)
    % order reduction
    [c, G, GI, E, id, dred] = nnHelper.reducePolyZono(c, G, GI, ...
        E, id, nrGen, obj.sensitivity);
    id_ = max(max(id), id_);

    % add approx error from reduction to error from approximation
    d = d + dred;

elseif max(obj.order) > 1 && ~evParams.sort_exponents
    % exponents remain equal with order = 1
    [E, G] = removeRedundantExponents(E, G);
end
end


function [G, GI, E, id, id_] = aux_addApproxError(obj, d, G, GI, E, id, id_, evParams)
% process approximation error

% save error bound for refinement
obj.refine_heu = d ...
    .* obj.do_refinement; % exclude neurons which should not be refined

D = diag(d);
D = D(:, d > 0);

if evParams.add_approx_error_to_GI
    GI = [GI, D];
else
    G = [G, D];
    E = blkdiag(E, eye(size(D, 2)));
    id = [id; (1:size(D, 2))' + id_];
    id_ = max(id);
end
end

% ------------------------------ END OF CODE ------------------------------
