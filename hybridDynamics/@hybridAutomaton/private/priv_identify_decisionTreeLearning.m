function node = priv_identify_decisionTreeLearning(data,options,varargin)
% priv_identify_decisionTreeLearning - This function constructs a state 
%   space partitioning via decision tree learning according to Section 3.2 
%   in [1]. The learning process is realized via a recursive function that 
%   repeateadly calls itself for growing the tree at a specific node.            
%
% Syntax:
%   node = priv_identify_decisionTreeLearning(data,options)
%   node = priv_identify_decisionTreeLearning(data,options,poly)
%
% Inputs:
%   -data:  struct storing the data points, the cluster assignments, and
%           other useful information that is required for decision tree 
%           learning
%   -options: a structure containing user-defined algorithm parameters
%
%       .minGini:  minimum Gini index at which the algorithm stops 
%                  growing the decision tree (see Section 3.2 in [1]).
%       .depth:    maximum depth of the subtree.
%       .noise:    array storing the noise estimate for each dimension
%       .ARX:      order for the AutoRegressive model with eXogenous
%                  inputs (ARX).
%
%   -poly:  polytope representing state space region of the current node 
%           (object of the polytope class from the CORA toolbox)      
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

    [poly,data_] = aux_initialization(data,options,nargin,varargin);

    % terminate decision tree growing if node is clean enough
    node.gini = aux_giniIndex(data.clus);
    [node.giniTraj,~,giniTraj] = aux_giniIndexTraj(data_.clus,data_.traj, ...
                                             data_.index,data_.x,options);
    jumps = aux_detectJumps(data,poly,options);

    if (node.gini < options.minGini || (node.giniTraj < options.minGini ...
        && giniTraj < options.minGini) || options.depth <= 0 || ...
        isempty(data_.x)) && isempty(jumps)

        node.data = data.index; node.poly = poly; 
        return;
    end

    % construct library of potential partitions (see Algorithm 1 in [1])
    if ~isempty(jumps)                      % jumps occur in the data
        
        % loop over all possible transitions
        fitted = []; surf = {}; points = {}; cnt = 1;

        for i = 1:max(data.clus)
            for j = 1:max(data.clus)
                [surf{cnt},fitted(cnt),points{cnt}] = ...
                            aux_identifyJumpingSurface(data,jumps,i,j,options);
                cnt = cnt + 1;
            end
        end

        % select surface which fits the most jumps within the noise
        [~,best] = max(fitted);
        node.type = 'jump';
        node.points = points{best};

    else                                    % no jumps occur in the data

        % loop over all possible transitions
        gini = []; fitted = []; surf = {}; type = {}; cnt = 1;

        for i = 1:max(data.clus)
            for j = 1:max(data.clus)
                if i ~= j
                    [surf{cnt},gini(cnt),fitted(cnt),type{cnt}] = ...
                               aux_identifyNonJumpingSurface(data_,i,j,options);
                    cnt = cnt + 1;
                end
            end
        end

        % catch the case where none of the identified surfaces is feasible
        if all(isinf(gini))
            node.data = data.index; node.poly = poly; return;
        end

        % select best surface based on Gini index and fitted changepoints
        best = aux_selectBestSurface(gini,fitted);

        node.type = type{best};
    end

    % construct child nodes
    node.surface = surf{best}; node.poly = poly; node.data = data.index;

    polyLeft = poly & polytope(node.surface.Ae,node.surface.be);
    polyRight = poly & polytope(-node.surface.Ae,-node.surface.be);
    
    tmp = node.surface.Ae * data.x - node.surface.be;
    ind = find(tmp <= 0); dataLeft = aux_subsetStruct(data,ind);
    ind = setdiff(1:length(tmp),ind); dataRight = aux_subsetStruct(data,ind);

    options.depth = options.depth - 1;

    if ~isempty(dataLeft.x) && ~isempty(dataRight.x)
        node.left = priv_identify_decisionTreeLearning(dataLeft,options,polyLeft);
        node.right = priv_identify_decisionTreeLearning(dataRight,options,polyRight);
    end
end


% Auxiliary functions -----------------------------------------------------

function [poly,data_] = aux_initialization(data,options,nargin,poly)
% intialize some auxiliary variables
   
     % initialize polytope
    if nargin <= 2
        poly = fullspace(size(data.x,1));
    else
        poly = normalizeConstraints(poly{1},'A');
    end

    % remove data points that are within the noise threshold from boundary
    if representsa(poly,'polytope')
       ind = find(all(abs(poly.A*data.x - poly.b) >= ...
                                            abs(poly.A*options.noise),1));
       data_ = aux_subsetStruct(data,ind);
    else
       data_ = data;
    end
end

function best = aux_selectBestSurface(gini,fitted)
% select the best switching surface based on the trajectory-based Gini 
% index and the number of changepoints fitted within the noise bound

    % sort surfaces according to the number of fitted changepoints
    [~,I] = sort(fitted,'descend'); 
    J1 = 1:numel(I);
    J1(I) = J1;

    % sort surfaces according to the trajectory-based Gini index
    [~,I] = sort(gini); 
    J2 = 1:numel(I);
    J2(I) = J2;

    % select the surface that performs best in both metrics
    J = J1 + J2;
    ind = find(~isinf(gini));
    [~,best] = min(J(ind));
    best = ind(best);
end

function gini = aux_giniIndex(clus)
% compute the Gini index according to Equation (14) in [1]

    points = zeros(max(clus));
    
    for j = 1:max(clus)
        points(j) = length(find(clus == j));
    end

    gini = 1-sum((points/sum(points)).^2);
end

function [gini,giniAll,giniTraj] = aux_giniIndexTraj(clus,traj,index,x,options)
% compute the Gini index under consideration of the trajectory information
% according to Equation (16) in [1]

    if isempty(clus)
        gini = inf; giniAll = [inf;inf]; giniTraj = inf; return;
    end

    % split data into trajectory segments
    seg = aux_trajectorySegments(traj,index);

    % compute Gini index for each segment
    cnt = 0; gini = 0; giniAll = [];

    for i = 1:length(seg)

        points = zeros(max(clus),1);
    
        for j = 1:max(clus)
            points(j) = length(find(clus(seg{i}) == j));
        end
    
        giniTmp = 1-sum((points/sum(points)).^2);
        giniAll = [giniAll,[giniTmp;sum(points)]];
        gini = gini + giniTmp*sum(points);
        cnt = cnt + sum(points);
    end

    gini = gini/cnt;

    % compute worst Gini index of all trajectory segments long enough
    if nargout > 2

        giniTraj = 0;

        for i = 1:length(seg)

            ind = seg{i};

            if length(ind) > 1

                len = sum(sqrt(sum((x(:,ind(2:end))-x(:,ind(1:end-1))).^2,1)));
    
                if len > 10*norm(options.noise)
                    giniTraj = max(giniTraj,giniAll(1,i));
                end
            end
        end
    end
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

function obj = aux_subsetStruct(obj,ind)
% return the subset of the struct specified by the given indices

   obj.x = obj.x(:,ind);
   obj.clus = obj.clus(ind);
   obj.index = obj.index(ind);
   obj.indTraj = obj.indTraj(ind);
   obj.traj = obj.traj(ind);
end

function jumps = aux_detectJumps(data,poly,options)
% detect indices where instant jumps occur in the data

    jumps = [];

    % loop over all trajectories segments
    seg = aux_trajectorySegments(data.traj,data.index);

    for i = 1:length(seg)

        % detect jumps based on the Euclidean distance
        x = data.x(:,seg{i});
        len = sqrt(sum((x(:,2:end) - x(:,1:end-1)).^2,1));
        ind = find(len > 2*norm(options.noise));

        % check if jumps are already covered by previous surfaces
        if representsa(poly,'polytope')
            ind = ind(any(abs(poly.A*x(:,ind) - poly.b) >= ...
                                        2*abs(poly.A*options.noise)));
        end

        jumps = [jumps,seg{i}(ind)];
    end
end

function [surf,fitted,points] = aux_identifyJumpingSurface(data,jumps,start,dest,options)
% identify a surface at which the trajectories jump

    % check if there exist any jumps for the current transition
    ind = find(data.clus(1:end-1) == start & data.clus(2:end) == dest);
    jumps = intersect(ind,jumps);

    if isempty(jumps)
        surf = []; points = []; fitted = -inf; return;
    end

    % determine points before the jumps
    jumps_ = [0;jumps]; xBefore = cell(length(jumps),1);

    for i = 2:length(jumps_)
        ind = jumps_(i-1)+1:jumps_(i);
        ind = intersect(ind,find(data.traj == data.traj(jumps_(i))));
        tmp = find(data.index(ind(2:end)) ~= data.index(ind(1:end-1)+1));
        if ~isempty(tmp)
            ind = ind(tmp(end):end);
        end
        xBefore{i-1} = data.x(:,ind);
    end

    % assemble data required for optimization
    xJump = data.x(:,jumps);

    dataOpt.x = xJump;
    dataOpt.xBefore = xBefore;
    dataOpt.noise = options.noise;
    
    % catch edge cases where there are not enough points for RANSAC alg.
    if size(xJump,2) < size(xJump,1)

        % determine jumping surface that fits all jumps perfectly
        dataDeg.x = horzcat(xBefore{:});
        dataDeg.clus = start*ones(size(dataDeg.x,2),1);
        dataDeg.noise = options.noise;

        Aeq = [xJump',-ones(size(xJump,2),1)];
        beq = zeros(size(Aeq,1),1);

        surf = aux_identifySwitchingSurface(dataDeg,start,0,Aeq,beq);

        w = [-surf.Ae';-surf.be];

    else

        % catch edge case where jump points form a perfect line
        [V,~,~] = svd(xJump-mean(xJump,2));
        w = V(:,end); w = [w;mean(w'*xJump)];

        xBefore = horzcat(xBefore{:});

        if sum(w(1:end-1)'*xBefore <= w(end)) <= size(xBefore,2)/2
            w = -w;
        end

        if any(abs(w(1:end-1)'*xJump - w(end)) > 2*abs(w(1:end-1)'*options.noise))
    
            % determine jumping surface using RANSAC optimization
            w = aux_ransacOptimization(xJump,dataOpt,@aux_costFunJump,options);
        end
    end

    % determine points before and after the fitted jumps
    [fitted,ind] = aux_costFunJump(w,dataOpt);
    
    fitted = -fitted;

    points.before = jumps(ind);
    points.after = jumps(ind)+1;

    surf = polytope([],[],w(1:end-1)',w(end));
end

function [surf,gini,fitted,type] = aux_identifyNonJumpingSurface(data,start,dest,options)
% identify a surface indicating a change in dynamics, which can either be a
% switching surface or a bouncing surface

    % check if there extis trasitions of this type in the data
    ind = find(data.clus(1:end-1) == start & data.clus(2:end) == dest,1);

    if isempty(ind)
        surf = []; gini = inf; type = 'none'; fitted = -inf; return;
    end

    % identify a swiching surface for this transition
    [surf1,cost1] = aux_identifySwitchingSurface(data,start,dest);

    % identify a bouncing surface for this transition
    [surf2,cost2] = aux_identifyBouncingSurface(data,start,dest,options);

    % select the surface with the lower cost
    if cost1 < cost2
        surf = surf1; type = 'switch';
    else
        surf = surf2; type = 'bounce';
    end

    % compute number of changepoints fitted within the noise bound
    ind = find(data.clus(1:end-1) == start & data.clus(2:end) == dest);
    x = data.x(:,ind);

    fitted = sum(abs(surf.Ae*x - surf.be) < abs(surf.Ae*options.noise));

    % compute the trajectory based Gini index
    tmp = surf.Ae * data.x - surf.be;
    indIn = find(tmp <= 0);
    indOut = setdiff(1:length(tmp),indIn);

    if isempty(indIn) || isempty(indOut)
        surf = []; gini = inf; type = 'none'; return;
    end

    [~,g1] = aux_giniIndexTraj(data.clus(indIn),data.traj(indIn),data.index(indIn));
    [~,g2] = aux_giniIndexTraj(data.clus(indOut),data.traj(indOut),data.index(indOut));

    gini = (sum(g1(1,:).*g1(2,:)) + sum(g2(1,:).*g2(2,:)))/(sum(g1(2,:))+sum(g2(2,:)));
end

function [surf,cost] = aux_identifySwitchingSurface(data,start,dest,varargin)
% identify a swiching surface where the dynamics changes and the trajectory
% crosses the surface

    % check if additional equality constraints are given
    Aeq = []; beq = [];

    if nargin > 3
        Aeq = varargin{1}; beq = varargin{2};
    end

    % determine data points for both sides of the surface
    xIn = data.x(:,data.clus == start);
    xOut = data.x(:,data.clus == dest);

    % compute an initial guess for the surface for faster optimization
    w0 = zeros(size(data.x,1)+1,1);

    if ~isempty(xIn) && ~isempty(xOut)
        mIn = median(xIn,2); mOut = median(xOut,2);
        w0 = mIn - mOut; w0 = w0/norm(w0);
        w0 = [w0;w0'*0.5*(mIn + mOut)];
    end

    % determine optimal switching surface via convex optimization
    options = optimoptions('fmincon','Display','off', ...
                            'SpecifyObjectiveGradient',true);

    w = fmincon(@(y) aux_costFunSetMembership(y,xIn,xOut), ...
                w0,[],[],Aeq,beq,[],[],[],options);

    % compute the cost (= number of misclassified points)
    cost = sum(w(1:end-1)'*xOut > w(end)) + sum(w(1:end-1)'*xIn < w(end));

    % convert surface to a constrained hyperplane object
    len = norm(w(1:end-1));
    surf = polytope([],[],w(1:end-1)'/len,w(end)/len);
end

function [surf,cost] = aux_identifyBouncingSurface(data,start,dest,options)
% identify a bouncing surface where the dynamics changes and the trajectory
% comes back into the same region

    xBorder = []; xDeg = []; clus = data.clus; indTraj = data.indTraj; 
    index = data.index; x = data.x; traj = data.traj;

    % split data into trajectory segments
    indInside = find(clus == start); indOutside = find(clus == dest);

    ind = find(indTraj(indInside(2:end)) ~= indTraj(indInside(1:end-1))+1 & ...
                    index(indInside(2:end)) ~= index(indInside(1:end-1))+1);
    ind = [0,ind,length(indInside)];

    segInside = cell(length(ind)-1,1);

    for i = 2:length(ind)
        segInside{i-1} = indInside(ind(i-1)+1:ind(i));
    end

    ind = find(indTraj(indOutside(2:end)) ~= indTraj(indOutside(1:end-1))+1 & ...
                    index(indOutside(2:end)) ~= index(indOutside(1:end-1))+1);
    ind = [0,ind,length(indOutside)];

    segOutside = cell(length(ind)-1,1);

    for i = 2:length(ind)
        segOutside{i-1} = indOutside(ind(i-1)+1:ind(i));
    end

    % detect segments that are connected
    data = {};
    indIn = 1:length(segInside);
    indOut = 1:length(segOutside);

    for i = 1:length(segInside)
        for j = 1:length(segOutside)
            in = segInside{i}; out = segOutside{j};
            if traj(in(1)) == traj(out(1)) && index(in(end)) + 1 == index(out(1))
                data{end+1}.x = [x(:,in),x(:,out)];
                data{end}.clus = [zeros(size(in));ones(size(out))];
                xBorder = [xBorder,x(:,in(end))];
                indIn = setdiff(indIn,i);
                indOut = setdiff(indOut,j);
            end
        end
    end

    % add data for segments that are not connected
    for i = 1:length(indIn)
        data{end+1}.x = x(:,segInside{indIn(i)});
        data{end}.clus = zeros(size(segInside{indIn(i)}));
    end

    for i = 1:length(indOut)
        data{end+1}.x = x(:,segOutside{indOut(i)});
        data{end}.clus = ones(size(segOutside{indOut(i)}));
        xDeg = [xDeg,data{end}.x];
    end

    % catch edge cases where there are not enough points for RANSAC alg.
    if isempty(xBorder) || (isempty(indIn) && isempty(indOut))
        surf = {}; cost = inf; return;
    elseif size(xBorder,2) < size(xBorder,1)
        dataDeg.x = [xDeg,x(:,clus == start)];
        dataDeg.clus = start*ones(size(dataDeg.x,2),1);
        Aeq = [xBorder',-ones(size(xBorder,2),1)];
        beq = zeros(size(Aeq,1),1);
        surf = aux_identifySwitchingSurface(dataDeg,start,dest,Aeq,beq);
        surf = polytope([],[],-surf.Ae,-(surf.be+1e-10));
        cost = aux_costFunBounce([surf.Ae';surf.be],data);
        return;
    end

    % catch edge case where switching points form a perfect line
    if options.ARX > 1
        [V,~,~] = svd(xBorder-mean(xBorder,2));
        c = V(:,end)'; d = mean(c*xBorder);

        if all(abs(c*xBorder - d) < abs(c*options.noise))
            if aux_costFunBounce([c';d],data) <= aux_costFunBounce(-[c';d],data)
                surf = polytope([],[],c,d); cost = 0; return;
            else
                surf = polytope([],[],-c,-d); cost = 0; return;
            end
        end
    end

    % identify bouncing surface via optimization
    [w,cost] = aux_ransacOptimization(xBorder,data,@aux_costFunBounce,options);
    
    surf = polytope([],[],w(1:end-1)',w(end));
end

function [w,cost] = aux_ransacOptimization(xBorder,data,costFun,options)
% determine best switching surface using the RANSAC algorithm

    N = 50;     % number of combinations that are tried

    n = size(xBorder,1);
    m = min(n,size(xBorder,2));
    cost = inf;

    % check if it is feasible to explore all combinations
    combs = nchoosek(size(xBorder,2),m);

    if combs < 10000
        combs = combinator(size(xBorder,2),m,'c');
        if size(combs,1) > N
            combs = combs(1:N,:);
        end
    else
        combs = zeros(N,m);
        for i = 1:N
            tmp = randperm(size(xBorder,2));
            combs(i,:) = tmp(1:m);
        end
    end

    % loop over all random picks
    for i = 1:size(combs,1)
        
        % randomly select n points to determine the swiching surface
        ind = combs(i,:);
        tmp = [ind,setdiff(1:size(xBorder,2),ind)];

        % learn switching surface using these points
        points = xBorder(:,ind);
        [V,~,~] = svd(points-mean(points,2));
        c = V(:,end);
        d = c'*xBorder(:,tmp(m:end));

        % determine shifts that are within the noise threshold
        tmp = d([false,abs(d(2:end) - d(1)) < 2*abs(c'*options.noise)]);

        if ~isempty(tmp)
            [~,ind] = sort(abs(tmp - d(1)),'ascend');
            d = [d(1),tmp(ind(min(length(ind),40)))];
        else
            d = d(1);
        end

        % determine optimal location of the surface by shifting it
        cost1 = costFun([c;d(1)],data);
        cost2 = costFun(-[c;d(1)],data);

        costAll = zeros(length(d),1);

        if cost1 < cost2    % use original direction of the surface
            costAll(1) = cost1;
            for j = 2:length(d)
                costAll(j) = costFun([c;d(j)],data);
            end
            [~,ind] = min(costAll);
            wTmp = [c;d(ind)];
            cTmp = costAll(ind);
        else                % use reverse direction of the surface
            costAll(1) = cost2;
            for j = 2:length(d)
                costAll(j) = costFun(-[c;d(j)],data);
            end
            [~,ind] = min(costAll);
            wTmp = -[c;d(ind)];
            cTmp = costAll(ind);
        end

        if cTmp < cost
            cost = cTmp; w = wTmp;
        end
    end
end

function [cost,fitted] = aux_costFunJump(w,data)
% cost for identifying jumping surfaces, where we give a reward for every
% jump that is fitted within the noise bound (see Equation (9) in [1])             

    fitted = [];

    % determine jumps that are fitted within the error bound
    ind = find(abs(w(1:end-1)'*data.x-w(end)) <= abs(w(1:end-1)'*data.noise));

    % check which jumps cross the hyperplane in the correct direction
    cost = 0;

    for i = 1:length(ind)
        in = sum(w(1:end-1)'*data.xBefore{ind(i)} <= w(end));
        out = size(data.xBefore{ind(i)},2) - in;
        if in > out
            cost = cost - 1;
            fitted = [fitted,ind(i)];
        end
    end
end

function cost = aux_costFunBounce(w,data)
% cost for bouncing off an surface, where the cost is the number of
% misclassified points (see Figure 6 in [1])

    c = w(1:end-1); d = w(end); cost = 0;

    % loop over all trajectories
    for i = 1:length(data)

        clus = data{i}.clus;

        % determine all crossings between trajectory and surface
        inside = c'*data{i}.x <= d;
        indIn = find(inside == 1);
        indOut = find(inside == 0);

        ind = find(inside(1:end-1) == 1 & inside(2:end) == 0);

        if isempty(ind)     % no correct crossing

            if length(unique(clus)) == 1
                cost = cost + length(find(clus(indOut) == 0));
            else
                cost = cost + length(find(clus(indIn) == 1));
                cost = cost + length(find(clus(indOut) == 0));
            end

        else                % at least one correct crossing

            % select the crossing with the lowest cost
            costTmp = zeros(length(ind),1);

            for j = 1:length(ind)

                before = 1:ind(j); 
                after = ind(j)+1:length(clus);

                costTmp(j) = costTmp(j) + length(find(clus(intersect(before,indIn)) == 1));
                costTmp(j) = costTmp(j) + length(find(clus(intersect(before,indOut)) == 0));

                costTmp(j) = costTmp(j) + length(find(clus(intersect(after,indIn)) == 0));
                costTmp(j) = costTmp(j) + length(find(clus(intersect(after,indOut)) == 1));
            end

            cost = cost + min(costTmp);
        end
    end
end

function [err,grad] = aux_costFunSetMembership(w,x1,x2)
% cost function for determining optimal switching surface by approximating
% the set-membership function with a sigmoid (see Equation (13) in [1])

    x1 = [x1;-ones(1,size(x1,2))]; x2 = [x2;-ones(1,size(x2,2))];

    s1 = 1./(1+exp(-w'*x1));
    err1 = sum((1-s1).^2);
    grad1 = sum(-2*(1-s1).*s1.*(1-s1).*x1,2);

    s2 = 1./(1+exp(-w'*x2));
    err2 = sum((0-s2).^2);
    grad2 = sum(-2*(0-s2).*s2.*(1-s2).*x2,2);

    err = err1 + err2;
    grad = grad1 + grad2;
end

% ------------------------------ END OF CODE ------------------------------
