function HA = priv_identify_decisionTree2automaton(tree,data,sys,noise)
% priv_identify_decisionTree2automaton - This function converts a state 
%   space partitioning represented by a decision tree to a hybrid automaton 
%   as described in Section 3.3 in [1].During the conversion, the 
%   information whether a partition corresponds to a switching surface, a 
%   bouncing surface, or a jumping surface is exploited.  
%
% Syntax:
%   HA = priv_identify_decisionTree2automaton(tree,data,sys,noise)
%
% Inputs:
%   -tree:    struct representing the decision tree
%   -data:    struct storing the data points, the cluster assignments, and
%             other useful information 
%   -sys:     cell-array storing the system dynamics for each cluster
%   -noise:   array storing the noise estimate for each dimension    
%
% Outputs:
%   -node:  struct storing the decision tree
%
% References:
%   [1] N. Kochdumper and et al. "Robust Identification of Hybrid Automata 
%       from Noisy Data", HSCC 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton/identify

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    n = size(data.x,1);

    % extract root nodes from the tree
    nodes = aux_extractRootNodes(tree);

    % catch the case with a single node
    if length(nodes) == 1
        tmp = unique(data.clus);
        [~,ind] = max(sum(data.clus' == tmp,2));
        HA = hybridAutomaton(location(nodes{1}.poly,transition(),sys{tmp(ind)}));
        return;
    end

    % remove redundancies in polytope representation
    nodes = aux_removeRedundanciesPolytope(nodes);

    % slightly bloat invariant regions for numerical stability
    tol = 1e-3;

    for i = 1:length(nodes)
        nodes{i}.poly = polytope(nodes{i}.poly.A,nodes{i}.poly.b+tol);
    end

    % determine which nodes correspond to overlapping regions
    modes = {};

    for i = 1:length(nodes)

        [dyn,seg,segDyn] = aux_determineDynamics(nodes{i},data,noise);
        
        % overlapping region
        if nodes{i}.giniTraj < nodes{i}.gini && length(dyn) > 1 && ...
                                    any(ismember(nodes{i}.type,'bounce'))

            nodes{i}.modes = [];

            for j = 1:length(dyn)
                modes{end+1,1} = struct('node',i,'dynamics',dyn(j));
                modes{end,1}.seg = segDyn{j};
                nodes{i}.modes = [nodes{i}.modes;length(modes)];
            end

        % non-overlapping region
        else
            [~,ind] = max(cellfun(@(x) length(x),segDyn));
            modes{end+1,1} = struct('node',i,'dynamics',dyn(ind));
            modes{end,1}.seg = seg;
            nodes{i}.modes = length(modes);
        end

        nodes{i}.seg = seg;
    end

    % determine potential target nodes for all outgoing surfaces
    for i = 1:length(nodes)

        nodes{i}.trans = cell(length(nodes{i}.poly.b),1);

        for j = 1:length(nodes)

            if i ~= j

                % loop over all potential outgoing transitions
                tmp = intersect(nodes{i}.surf,nodes{j}.surf);
    
                for k = 1:length(tmp)
    
                    % construct guard set
                    ind1 = find(nodes{i}.surf == tmp(k));
                    A = nodes{i}.poly.A; b = nodes{i}.poly.b;
                    c = A(ind1,:); d = b(ind1);
                    A(ind1,:) = []; b(ind1,:) = [];
    
                    ind2 = find(nodes{j}.surf == tmp(k));
                    A_ = nodes{j}.poly.A; b_ = nodes{j}.poly.b;
                    c_ = A_(ind2,:);
                    A_(ind2,:) = []; b_(ind2,:) = [];
                    A = [A;A_]; b = [b;b_];
    
                    guard = polytope(A,b,c,d);

                    % add transition to list
                    if ~representsa(guard,'emptySet') && norm(c + c_) < eps
                        tran = struct('guard',guard,'target',j);
                        nodes{i}.trans{ind1} = [nodes{i}.trans{ind1};tran];
                    end
                end
            end
        end
    end

    % remove overlapping guard sets
    nodes = aux_removeOverlappingGuards(nodes);

    % construct transitions for each mode
    for i = 1:length(modes)
        
        node = nodes{modes{i}.node};
        modes{i}.trans = []; modes{i}.transInd = []; modes{i}.transPoi = [];

        % loop over all surfaces of the polytope
        for j = 1:length(node.poly.b)

            % construct guard set for this surface
            ind = setdiff(1:length(node.poly.b),j);

            guard = polytope(node.poly.A(ind,:),node.poly.b(ind), ...
                             node.poly.A(j,:),node.poly.b(j)-tol);

            % distinguish between the different types of surfaces
            if strcmp(node.type{j},'jump') && strcmp(node.side{j},'left')
                
                % identify reset function z = A*x + c (Equation (17) in [1])
                [A,c] = aux_identifyResetFunction(node.points{j},data);

                % check which nodes are reachable by the jump
                for k = 1:length(nodes)
                    
                    % construct guard set for jumping into current node
                    % (see Proposition 1 in [1])
                    C = nodes{k}.poly.A*A;
                    d = nodes{k}.poly.b - nodes{k}.poly.A*c;

                    C = [guard.A;C]; d = [guard.b;d];

                    guardJump = polytope(C,d,guard.Ae,guard.be+tol);

                    % check if it is possible to jump into the current node
                    if ~representsa(guardJump,'emptySet')
                        
                        % determine the potential target modes for this transition
                        targets = nodes{k}.modes;

                        % select the most likely target mode if there are multiple
                        dyn = cellfun(@(x) x.dynamics,modes(targets));
                        ind = aux_dynamicsAfterJump(node.points{j}.after,nodes{k},data,dyn);
                        targets = targets(ind);
        
                        % construct transition object
                        reset = linearReset(A,[],c);
                        transTemp = transition(guardJump,reset,targets);
        
                        modes{i}.trans = [modes{i}.trans,transTemp];

                        modes{i}.transInd = [modes{i}.transInd;0];
                        modes{i}.transPoi = [modes{i}.transPoi;0];
                    end
                end

            elseif strcmp(node.type{j},'bounce') && strcmp(node.side{j},'left') ...
                                                            && length(node.modes) > 1

                % determine the potential target modes for this transition
                targets = setdiff(node.modes,i);

                % select the most likely target mode if there are multiple
                dyn = cellfun(@(x) x.dynamics,modes(targets));
                [ind,success,points] = aux_dynamicsAfterBouncing(modes{i},node,node.poly,j,data,dyn);
                targets = targets(ind);

                % construct transition object
                if success
                    reset = linearReset(eye(n));
                    transTemp = transition(guard,reset,targets);
    
                    modes{i}.trans = [modes{i}.trans,transTemp];

                    modes{i}.transInd = [modes{i}.transInd;j];
                    modes{i}.transPoi = [modes{i}.transPoi;points];
                end
            end

            % create default switching surfaces by looping over all guard 
            % sets for this surface
            if ~(strcmp(node.type{j},'jump') && strcmp(node.side{j},'left'))
                for k = 1:length(node.trans{j})
    
                    % determine the potential target modes for this transition
                    targets = nodes{node.trans{j}(k).target}.modes;
    
                    % select the most likely target mode if there are mulitple
                    if length(targets) > 1
                        dyn = cellfun(@(x) x.dynamics,modes(targets));
                        ind = aux_dynamicsAfterCrossing(node, ...
                                        nodes{node.trans{j}(k).target},data,dyn);
                        targets = targets(ind);
                    end
    
                    % construct transition object
                    reset = linearReset(eye(n));
                    transTemp = transition(node.trans{j}(k).guard,reset,targets);
        
                    modes{i}.trans = [modes{i}.trans,transTemp];

                    modes{i}.transInd = [modes{i}.transInd;0];
                    modes{i}.transPoi = [modes{i}.transPoi;0];
                end
            end
        end
    end

    % remove conflicting bouncing surfaces
    modes = aux_removeConflictingSurfaces(modes);

    % construct hybrid automaton (with relearned dynamics)
    loc = [];

    for i = 1:length(modes)
        node = nodes{modes{i}.node};
        dyn = aux_identifyLinearSys(data,modes{i}.seg, ...
                                sys{modes{i}.dynamics},noise);
        loc = [loc; location(node.poly,modes{i}.trans,dyn)];
    end

    HA1 = hybridAutomaton(loc);

    % construct hybrid automaton (without relearned dynamics)
    loc = [];

    for i = 1:length(modes)
        node = nodes{modes{i}.node};
        dyn = sys{modes{i}.dynamics};
        loc = [loc; location(node.poly,modes{i}.trans,dyn)];
    end

    HA2 = hybridAutomaton(loc);

    % select best hybrid automaton
    err1 = aux_predictionErrorAutomaton(HA1,data);
    err2 = aux_predictionErrorAutomaton(HA2,data);

    if err1 < err2
        HA = HA1;
    else
        HA = HA2;
    end
end


% Auxiliary functions -----------------------------------------------------

function [rootNodes,cnt] = aux_extractRootNodes(tree,varargin)
% extract the root nodes from a decision tree

    % initialization
    if nargin == 1
        surf.id = [];
        surf.type = {};
        surf.side = {};
        surf.points = {};
        surf.cnt = 0;
    else
        surf = varargin{1};
    end

    % check if the current nodes has children
    if isfield(tree,'left') && isfield(tree,'right')

        surf.cnt = surf.cnt + 1;
        surf.id = [surf.id;surf.cnt];
        surf.type = [surf.type,{tree.type}];
        surf.side = [surf.side,{'left'}];
        surf.points = [surf.points,{[]}];
        if isfield(tree,'points')
            surf.points{end} = tree.points;
        end

        [nodes1,cnt] = aux_extractRootNodes(tree.left,surf);

        surf.cnt = cnt;
        surf.side{end} = 'right';

        [nodes2,cnt] = aux_extractRootNodes(tree.right,surf);

        rootNodes = [nodes1,nodes2];

    % if the nodes does not have children, its a root node
    else
        cnt = surf.cnt;
        tree.surf = surf.id;
        tree.type = surf.type;
        tree.side = surf.side;
        tree.points = surf.points;
        rootNodes = {tree};
    end
end

function [dyn,segAll,segDyn] = aux_determineDynamics(node,data,noise)
% determine which dynamics occur in the current mode

    % split the data into single trajectory segments
    clus = data.clus(node.data); traj = data.traj(node.data);
    index = data.index(node.data); x = data.x(:,node.data);

    seg = aux_trajectorySegments(traj,index);

    % determine dynamics for each segment via majority vote
    dyn = zeros(size(seg)); small = zeros(size(seg)); segAll = seg;
    pointsAll = zeros(max(clus),1);

    for i = 1:length(seg)

        segAll{i} = index(seg{i});

        % check if the segment is smaller than the noise
        if all(all(abs(x(:,seg{i}) - mean(x(:,seg{i}),2)) <= noise))
            small(i) = 1;
        end

        % perform a majority vote to obtain the dynamics
        clus_ = clus(seg{i});
        points = zeros(max(clus),1);

        for j = 1:length(points)
            points(j) = sum(clus_ == j);
        end

        [~,dyn(i)] = max(points);

        pointsAll = pointsAll + points;
    end
    
    % disregard segments that are smaller than the noise
    ind = find(small == 1);

    if length(ind) == length(small)
        [~,dyn] = max(pointsAll); segDyn = {segAll}; return;
    end

    ind = setdiff(1:length(small),ind);
    dyn = dyn(ind); segAll = segAll(ind);

    % determine all dynamics that occur in the segments
    [dyn,~,ind] = unique(dyn);

    segDyn = cell(length(dyn),1);

    for i = 1:length(dyn)
        segDyn{i} = segAll(ind == i);
    end
end

function [A,c] = aux_identifyResetFunction(points,data)
% identify the paramters for the reset function z = A*x + c according to 
% Equation (17) in [1]

    % collect data
    X1 = data.x(:,points.before);
    X1 = [X1;ones(1,size(X1,2))];
    X2 = data.x(:,points.after);

    % apply leave-one-out cross-validation to determine optimal rank
    if size(X1,2) > size(X1,1)

        err = zeros(size(X1,1),1);
        ind = unique(randi(size(X1,2),[10,1]));
    
        for i = 1:length(ind)
            index = setdiff(1:size(X1,2),ind(i));
            [V,S,W] = svd(X1(:,index),'econ');
            for j = 1:sum(diag(S > 0))
                V_ = V(:,1:j); S_ = S(1:j,1:j); W_ = W(:,1:j);
                A = X2(:,index)*W_*diag(1./diag(S_))*V_';
                err(j) = err(j) + norm(X2(:,ind(i)) - A*X1(:,ind(i)));
            end
        end
    
        [~,rank] = min(err);
    else
        rank = size(X1,1);
    end

    % determine best fit via Dynamic Mode Decomposition
    [V,S,W] = svd(X1,'econ');

    rank = min(rank,sum(diag(S) > 0));
    V = V(:,1:rank); S = S(1:rank,1:rank); W = W(:,1:rank);

    A = X2*W*diag(1./diag(S))*V';

    % extract optimal parameters
    c = A(:,end); A = A(:,1:end-1);
end

function [target,success] = aux_dynamicsAfterJump(jump,node,data,dyn)
% determine the most likely dynamics after jumping to another mode

    points = zeros(size(dyn));
    
    % loop over all trajectory segments in the new mode
    for i = 1:length(node.seg)

        % check if the segments is reached after the jump
        if ismember(node.seg{i}(1),jump)

            % compute number of correctly classified points
            for k = 1:length(dyn)
                points(k) = sum(data.clus(node.seg{i}) == dyn(k));
            end
        end
    end

    % select most likely dynamics via majority vote
    [~,target] = max(points);

    % try to find trajectories corresponding to jumps if no point was found
    success = any(points > 0);

    if ~success

        % loop over all jumps
        for j = 1:length(jump)

            cost = inf*ones(length(node.seg),1);

            % loop over all segments
            for i = 1:length(node.seg)

                if data.traj(jump(j)) == data.traj(node.seg{i}(1)) && ...
                        node.seg{i}(1) > jump(j)
                    cost(i) = node.seg{i}(1) - jump(j);
                end
            end

            % determine best segment
            [~,ind] = min(cost);

            % compute number of correctly classified points
            if ~isinf(cost(ind))
                for k = 1:length(dyn)
                    points(k) = sum(data.clus(node.seg{ind}) == dyn(k));
                end 
            end
        end

        % select most likely dynamics via majority vote
        [~,target] = max(points);

        success = any(points > 0);
    end
end

function target = aux_dynamicsAfterCrossing(node1,node2,data,dyn)
% determine the most likely dynamics after crossing to another mode
    
    points = zeros(size(dyn));

    % loop over all segment combination
    for i = 1:length(node1.seg)
        for j = 1:length(node2.seg)
            
            % check if the segments belong together
            if data.traj(node1.seg{i}(1)) == data.traj(node2.seg{j}(1)) && ...
                node1.seg{i}(end) < length(data.index) && ...
                data.index(node1.seg{i}(end)+1) == data.index(node2.seg{j}(1))

                % compute number of correctly classified points
                for k = 1:length(dyn)
                    points(k) = sum(data.clus(node1.seg{i}) == dyn(k));
                end
            end
        end
    end

    % select most likely dynamics via majority vote
    [~,target] = max(points);
end

function [target,success,score] = aux_dynamicsAfterBouncing(mode,node,poly,guard,data,dyn)
% determing the most likely dynamics after bouncing of a surface

    points = zeros(size(dyn));

    % determine segments in the current node that bounce of the surface
    index = [];

    for i = 1:length(node.seg)

        % check if the segment bounces off the current surface
        tmp = abs(poly.A * data.x(:,node.seg{i}(1)) - poly.b);
        [~,ind] = min(tmp);

        if ind == guard && data.indTraj(node.seg{i}(1)) > 1
            index = [index;i];
        end
    end
    
    seg = node.seg(index);

    % loop over all segments in the current mode
    for i = 1:length(mode.seg)

        % check if the segment bounces off the current surface
        tmp = abs(poly.A * data.x(:,mode.seg{i}(end)) - poly.b);
        [~,ind] = min(tmp);

        if ind == guard && mode.seg{i}(end) < length(data.traj) && ...
                data.traj(mode.seg{i}(end)+1) == data.traj(mode.seg{i}(end))

            % loop over all segments in the current node
            for j = 1:length(seg)
                
                % check if segments can be connected
                if data.traj(mode.seg{i}(1)) == data.traj(seg{j}(1)) && ...
                   data.indTraj(seg{j}(1)) > data.indTraj(mode.seg{i}(1))

                    outside = data.indTraj(seg{j}(1)) - ...
                                        data.indTraj(mode.seg{i}(end));
                    
                    for k = 1:length(dyn)
                        inside = sum(data.clus(seg{j}) == dyn(k));
                        if inside > outside
                            points(k) = points(k) + inside - outside;
                        end
                    end
                end
            end
        end
    end

    % determine most likely dynamcis via majority vote
    [score,target] = max(points);
    
    success = any(points > 0);
end

function nodes = aux_removeOverlappingGuards(nodes)
% remove overlapping guard sets to avoid non-determinism and improve
% numerical stability of the resulting hybrid automaton

    % loop over all nodes and surfaces of the corresponding polytope
    for i = 1:length(nodes)
        for j = 1:length(nodes{i}.trans)

            remove = []; trans = nodes{i}.trans{j};
  
            % loop over all combinations of guard sets for this surface
            for k = 1:length(trans)
                for h = 1:length(trans)
                    if k ~= h

                        % check if the guard set is already covered by
                        % another guard set
                        poly = polytope(trans(k).guard.A,trans(k).guard.b);

                        if contains(poly,trans(h).guard,'exact',1e-6)
                            remove = [remove;h];
                        end
                    end
                end
            end

            % remove all guard sets that are already covered by others
            if ~isempty(remove)
                keep = setdiff(1:length(trans),unique(remove));
                nodes{i}.trans{j} = nodes{i}.trans{j}(keep);
            end
        end
    end
end

function nodes = aux_removeRedundanciesPolytope(nodes)
% remove redundant halfspaces from the polytopes that represent the
% invariant regions 

    % loop over all root nodes
    for i = 1:length(nodes)

        tmp = nodes{i}.poly;

        % remove redundant halfspaces
        poly = compact(nodes{i}.poly,'aligned');

        % check which halfspaces were removed
        matOrig = [tmp.A,tmp.b]; matNew = [poly.A,poly.b];
        ind = zeros(size(matNew,1),1);

        for j = 1:size(matNew,1)
            [~,ind(j)] = min(sum(abs(matOrig - matNew(j,:)),2));
        end

        % remove stored informations for the removed halfspace 
        nodes{i}.surf = nodes{i}.surf(ind);
        nodes{i}.type = nodes{i}.type(ind);
        nodes{i}.side = nodes{i}.side(ind);
        nodes{i}.points = nodes{i}.points(ind);

        nodes{i}.poly = polytope(tmp.A(ind,:),tmp.b(ind));
    end
end

function modes = aux_removeConflictingSurfaces(modes)
% remove bouncing surfaces that are conflicting because after bouncing from
% one surface puts the state directly on the other one -> endless loop
% wihtout any progress in time

    % initialize list with indices of surfaces that are removed
    for i = 1:length(modes)
        modes{i}.remove = [];
    end

    % loop over all mode combinations
    for i = 1:length(modes)
        for j = i+1:length(modes)

            % check if modes are generated from the same decision tre node
            if modes{i}.node == modes{j}.node

                % determine common bouncing surfaces
                trans1 = modes{i}.transInd(modes{i}.transInd > 0);
                trans2 = modes{j}.transInd(modes{j}.transInd > 0);
                trans = intersect(trans1,trans2);

                % loop over all common bouncing surfaces
                for k = 1:length(trans)

                    ind1 = find(modes{i}.transInd == trans(k));
                    ind2 = find(modes{j}.transInd == trans(k));

                    % check if the bouncing surfaces are conflicting
                    if modes{i}.trans(ind1).target == j || ...
                            modes{j}.trans(ind2).target == i

                        % remove surface that fits less points
                        if modes{i}.transPoi(ind1) < modes{j}.transPoi(ind2)
                            modes{i}.remove = [modes{i}.remove,ind1];
                        else
                            modes{j}.remove = [modes{j}.remove,ind2];
                        end
                    end
                end
            end
        end
    end

    % remove all conflicting surfaces that have beed determined
    for i = 1:length(modes)
        if ~isempty(modes{i}.remove)
            ind = unique(modes{i}.remove);
            ind = setdiff(1:length(modes{i}.trans),ind);
            modes{i}.trans = modes{i}.trans(ind);
        end
    end
end

function A = aux_dynamicModeDecomposition(X1,X2,U)
% compute the matrix X2 = A*X1 that best fits the data

    % extend matrices by entries for constant inputs and inputs
    X1 = [X1; ones(1,size(X1,2))];

    if ~isempty(U)
        X1 = [X1; U];
    end

    % singular value decomposition
    [V,S,W] = svd(X1,'econ');
    
    % construct matrices with different rank
    rankMax = sum(diag(S) > 0);
    A = cell(length(diag(S)),1);

    for rank = 1:rankMax

        % reduce rank by removing the smallest singular values
        if ~isempty(rank) && rank < size(S,1)
            V_ = V(:,1:rank); S_ = S(1:rank,1:rank); W_ = W(:,1:rank);
        else
            V_ = V; S_ = S; W_ = W;
        end
    
        % compute resulting system matrix A
        A{rank} = X2*W_*diag(1./diag(S_))*V_';
    end
end

function sys = aux_identifyLinearSys(data,seg,default,noise)
% identify a linear system best fitting the data

    % split segments at jumps
    jumps = true;

    while jumps

        jumps = false;
    
        for i = 1:length(seg)
            len = sqrt(sum((data.x(:,seg{i}(2:end))-data.x(:,seg{i}(1:end-1))).^2,1));
            ind = find(len > 2*norm(noise));
            if ~isempty(ind)
                ind = [ind,length(seg{i})];
                for j = 2:length(ind)
                    seg{end+1} = seg{i}(ind(j-1)+1:ind(j));
                end
                seg{i} = seg{i}(1:ind(1));
                jumps = true;
                break;
            end
        end
    end

    % bring data to the correct format
    data_ = cell(length(seg),1);
    x = []; dx = []; u = [];

    for i = 1:length(seg)

        if length(seg{i}) > 1

            % split data into trajectory segments
            data_{i}.x = data.x(:,seg{i});
            data_{i}.t = data.t(seg{i});
    
            data_{i}.u = [];
    
            if isfield(data,'u')
                data_{i}.u = data.u(:,seg{i});
            end

            % recompute derivative to make sure only points from the
            % current segment are used
            data_{i}.dx = aux_derivative(data_{i}.t,data_{i}.x')';
    
            % add data from current segment to overall data
            x = [x,data_{i}.x]; dx = [dx,data_{i}.dx]; u = [u,data_{i}.u];
        end
    end

    if isempty(x)
        sys = default; return;
    end

    data = data_;
    data = data(~cellfun('isempty',data));

    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(x,dx,u);

    % select the matrix that best fits the data
    errBest = inf;

    for i = 1:length(Alist)
        err = 0;
        if ~isempty(Alist{i})
            for j = 1:length(data)
                if length(data{j}.t) > 1
                    errTemp = aux_computeError(Alist{i},data{j});
                    err = err + errTemp;
                end
            end
            if err < errBest
                Aall = Alist{i};
                errBest = err;
            end
        end
    end 

    % construct linear system object
    n = size(Aall,1);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 
    
    sys = linearSys(A,B,c);
end

function err = aux_computeError(Aall,data)
% compute the error between the system approximation and the real data

    % construct linear system object
    n = size(Aall,1);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 
    dt = data.t(2) - data.t(1);

    % loop over all time steps and integrate using Runge-Kutta-4
    x = zeros(size(data.x));
    x(:,1) = data.x(:,1);

    if isempty(data.u)
        for i = 1:size(x,2)-1
            k1 = A*x(:,i) + c;
            k2 = A*(x(:,i) + dt*k1/2) + c;
            k3 = A*(x(:,i) + dt*k2/2) + c;
            k4 = A*(x(:,i) + dt*k3) + c;
            x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    else
        for i = 1:size(x,2)-1
            temp = B*data.u(:,i) + c;

            k1 = A*x(:,i) + temp;
            k2 = A*(x(:,i) + dt*k1/2) + temp;
            k3 = A*(x(:,i) + dt*k2/2) + temp;
            k4 = A*(x(:,i) + dt*k3) + temp;
            x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    end

    % compute error
    err = sum(sqrt(sum((x-data.x).^2,1)));
end

function seg = aux_trajectorySegments(traj,index)
% split the data into continuous trajectory segments

    seg = {};

    ind1 = find(traj(1:end-1) ~= traj(2:end));
    ind1 = [0,ind1,length(traj)];

    for i = 1:length(ind1)-1

        tmp1 = ind1(i)+1:ind1(i+1);
        ind2 = find(index(tmp1(2:end)) ~= index(tmp1(1:end-1))+1);
        ind2 = [0,ind2,length(tmp1)];

        for j = 1:length(ind2)-1
            seg{end+1} = tmp1(ind2(j)+1:ind2(j+1));
        end
    end
end
function dy = aux_derivative(x,y)
    
    % number of subintervals
    N = length(x)-1;
    
    % preallocates vector to store derivative
    dy = zeros(size(y));
    
    % approximates derivative at lower bound using forward difference
    dy(1,:) = (y(2,:)-y(1,:))/(x(2)-x(1));
    
    % approximates derivative at upper bound using backward difference
    dy(N+1,:) = (y(N+1,:)-y(N,:))/(x(N+1)-x(N));
    
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(i,:) = (y(i+1,:)-y(i-1,:))/(x(i+1)-x(i-1));
    end
end

function err = aux_predictionErrorAutomaton(HA,data)
% compute the error for the current hybrid automaton

    % convert data to a trajectory based structure
    traj = cell(max(data.traj),1);

    for i = 1:length(traj)
        ind = find(data.traj == i);
        traj{i}.x = data.x(:,ind);
        traj{i}.t = data.t(ind);
        if isfield(data,'u')
            traj{i}.u = data.u(:,ind);
        end
    end

    % initialization
    err = 0;
    dt = traj{1}.t(2) - traj{1}.t(1);
    inputs = false;
    
    if isfield(traj{1},'u')
        inputs = true;
    end

    % loop over all trajectories
    for j  = 1:length(traj)

        % determine potential initial modes 
        initMode = [];

        for i = 1:length(HA.location)
            if contains(HA.location(i).invariant,traj{j}.x(:,1))
                initMode = [initMode;i];
            end
        end

        % loop over all potential initial modes
        errTmp = zeros(length(initMode),1);

        for k = 1:length(initMode)

            loc = HA.location(initMode(k));
            sys = loc.contDynamics; trans = loc.transition;
            A_ = sys.A; B_ = sys.B; c_ = sys.c; 

            % loop over all time steps and integrate using Runge-Kutta-4
            x = zeros(size(traj{j}.x));
            x(:,1) = traj{j}.x(:,1);
        
            if ~inputs
                for i = 1:size(x,2)-1

                    % select correct dynamics for the current point
                    if i > 1
                        for t = 1:length(trans)
                            if aux_guardCrossed(trans(t).guard,x(:,i),x(:,i-1))
                                x(:,i) = trans(t).reset.A*x(:,i) + trans(t).reset.c;
                                loc = HA.location(trans(t).target);
                                sys = loc.contDynamics; trans = loc.transition;
                                A_ = sys.A; c_ = sys.c; break; 
                            end
                        end
                    end

                    % integrate using Runge-Kutta-4
                    k1 = A_*x(:,i) + c_;
                    k2 = A_*(x(:,i) + dt*k1/2) + c_;
                    k3 = A_*(x(:,i) + dt*k2/2) + c_;
                    k4 = A_*(x(:,i) + dt*k3) + c_;

                    x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
                end
            else
                for i = 1:size(x,2)-1

                    % select correct dynamics for the current point
                    if i > 1
                        for t = 1:length(trans)
                            if aux_guardCrossed(trans(t).guard,x(:,i),x(:,i-1))
                                x(:,i) = trans(t).reset.A*x(:,i) + trans(t).reset.c;
                                loc = HA.location(trans(t).target);
                                sys = loc.contDynamics; trans = loc.transition;
                                A_ = sys.A; B_ = sys.B; c_ = sys.c; break;
                            end
                        end
                    end
                    
                    % integrate using Runge-Kutta-4
                    temp = B_*traj{j}.u(:,i) + c_;

                    k1 = A_*x(:,i) + temp;
                    k2 = A_*(x(:,i) + dt*k1/2) + temp;
                    k3 = A_*(x(:,i) + dt*k2/2) + temp;
                    k4 = A_*(x(:,i) + dt*k3) + temp;
                    
                    x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
                end
            end
    
            % compute error
            errTmp(k) = sum(sqrt(sum((x-traj{j}.x).^2,1)));
        end

        % select minimum error from all potential initial modes
        err = err + min(errTmp);
    end
end

function res = aux_guardCrossed(guard,x1,x2)
% check if the trajectory crossed a guard set represented as a constrained
% hyperplane

    if guard.Ae*x1 <= guard.be
        res = false; return;
    elseif isempty(guard.A)
        res = true; return;
    else

        % comptue point on the hyperplane
        d = (guard.be - guard.Ae*x1)/(guard.Ae*x2 - guard.Ae*x1);
        p = x1 + d*(x2 - x1);

        % check if the point satisfies the constraints
        res = all(guard.A * p <= guard.b);
    end
end

% ------------------------------ END OF CODE ------------------------------
