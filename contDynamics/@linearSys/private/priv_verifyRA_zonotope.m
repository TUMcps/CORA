function [res,R] = priv_verifyRA_zonotope(linsys,params,options,spec)
% priv_verifyRA_zonotope - verification of reach-avoid specifications for linear
%    systems using [1, Alg. 3]
%
% Syntax:
%    [res,R] = priv_verifyRA_zonotope(linsys,params,options,spec)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - settings for reachability analysis
%    spec - object of class specification
%
% Outputs:
%    res - true/false (true if specifications verified, otherwise false)
%    R - reachSet object, outer-/inner-approximative depending on res
%
% References:
%    [1] M. Wetzlinger et al. "Fully automated verification of linear
%        systems using inner-and outer-approximations of reachable sets",
%        TAC, 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/verify

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       16-August-2021
% Last update:   08-November-2022 (MW, proper integration in master)
%                29-September-2024 (TL, init R if desired)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % extract safe sets and unsafe sets from specifications
    [params.safeSet,params.unsafeSet] = aux_getSetsFromSpec(spec);

    % input argument pre-processing
    [params,options] = validateOptions(linsys,params,options);
    
    % estimate the required max. allowed error with simulations
    emax = aux_initError(linsys,params);

    % verification algorithm using adaptive reachability algorithm
    options.verify = true;
    options.linAlg = 'adaptive';

    % init reachable set containers
    if nargout == 2
        RtimePoint.set = {}; RtimePoint.time = {}; RtimePoint.error = [];
        RtimeInt.set = {}; RtimeInt.time = {}; RtimeInt.error = [];
    end
    
    % refinement loop
    cnt = 0;
    while true
        
        % increment counter
        cnt = cnt + 1;
       
        % init/reset distance values
        d_Go = -Inf; d_Gi = -Inf; d_Fo = -Inf; d_Fi = Inf;
        ind_Go = 0; ind_Gi = 0; ind_Fo = 0; ind_Fi = 0;
        
        % compute outer-approximation of the reachable set
        options.error = emax;
        if options.verbose
            fprintf("Iteration %i, current error: %d\n", cnt, options.error);
        end
        [timeInt,timePoint,res,savedata] = priv_reach_adaptive(linsys,params,options);

        % rewrite information about specification satisfaction into params
        for i=1:numel(params.safeSet)
            params.safeSet{i}.time_ = savedata.safeSet{i}.time_;
        end
        for i=1:numel(params.unsafeSet)
            params.unsafeSet{i}.time_ = savedata.unsafeSet{i}.time_;
        end
        savedata = rmiffield(savedata,{'safeSet','unsafeSet'});
        options.savedata = savedata;
        
        % save first time-point solution (same in every iteration)
        if cnt == 1
            RtimePoint.set = [RtimePoint.set; timePoint.set(1)];
            RtimePoint.time = [RtimePoint.time; {params.tStart}];
            RtimePoint.error = [RtimePoint.error; timePoint.error(1)];
        end
        
        % already falsified or (completely!) verified in priv_reach_adaptive
        allSpecTimes = [cellfun(@(S) S.time_,params.safeSet,'UniformOutput',false); ...
                        cellfun(@(S) S.time_,params.unsafeSet,'UniformOutput',false)];
        if ~res || all(representsa(vertcat(allSpecTimes{:}),'emptySet',0))
            % create reachSet object if desired
            if nargout == 2
                R = aux_initReachableSet(RtimePoint,RtimeInt);
            end
            return
        end
        
        % generate sequence of time intervals
        time = repmat(cell2mat(timePoint.time(2:end)),1,2);
        time(:,1) = [0; time(1:end-1,1)];
        % flag to see which sets can be saved in R
        allSAT = false(numel(timeInt.set),numel(params.safeSet)+numel(params.unsafeSet));
        allSATbefore = true(numel(timeInt.set),1);
        
        
        % clear variables for inner-approximation between runs
        R_i = {}; ind_Ri = [];
        
        % check if the reachable set is inside the safe sets       
        for i = 1:numel(params.safeSet)
            for j = 1:numel(timeInt.set)

                if isIntersecting(params.safeSet{i}.time_,time(j,:))
                    % outer-approximation
                    d = aux_containmentCheckZono(timeInt.set{j},params.safeSet{i}.set);
                    
                    if d < 0
                        % outer-approximation is fully contained in the
                        % safe set -> extend already verified time intervals
                        params.safeSet{i}.time_ = ...
                            setdiff(params.safeSet{i}.time_,time(j,:));

                        % append to reachable set
                        allSAT(j,i) = true;
                        allSATbefore(j) = false;
                    end
                    
                    if d > d_Go
                        d_Go = d; ind_Go = j; 
                    end
                    
                    % quick check for inner-approximation
                    if params.safeSet{i}.fastInner && d > timeInt.error(j)
                        % the distance d from the over-approximation to
                        % the safe set is larger than the contained
                        % error so that the inner-approximation will
                        % never be fully contained in the safe set
                        % -> provably unsafe in this special case
                        res = false; 
                        if nargout == 2
                            R = aux_initReachableSet(RtimePoint,RtimeInt);
                        end
                        return;
                    end
                    
                    % inner-approximation
                    if d > timeInt.error(j)
                        % outer-approximation not fully contained in the
                        % safe set -> check whether inner-approximation
                        % intersects the safe set
                        [R_i,ind_Ri,ind_now] = aux_compInnerApprox(...
                            R_i,ind_Ri,timeInt.set,timeInt.error,j);
                        d = aux_containmentCheckConZono(R_i{ind_now},params.safeSet{i}.set);
                        if d > 0
                            % inner-approximation not fully contained in
                            % the safe set -> provably unsafe
                            res = false; 
                            if nargout == 2
                                R = aux_initReachableSet(RtimePoint,RtimeInt);
                            end
                            return
                        end
                    else
                        d = d - timeInt.error(j);
                    end

                    if d > d_Gi
                    	d_Gi = d; ind_Gi = j; 
                    end
                else
                    % append to reachable set
                    allSAT(j,i) = true;
                end
            end
        end
        
        % check if the reachable set is inside the unsafe sets
        for i = 1:numel(params.unsafeSet)
            for j = 1:numel(timeInt.set)
                
                % check if specification is (still) active
                if isIntersecting(params.unsafeSet{i}.time_,time(j,:))
                
                    % fast intersecting check for speed-up
                    if aux_intersectionCheckFast(timeInt.set{j},params.unsafeSet{i})

                        % outer-approximation
                        d = aux_distanceIntersectionZono(timeInt.set{j},params.unsafeSet{i}.set);

                        if d > d_Fo
                            d_Fo = d; ind_Fo = j; 
                        end

                        if params.unsafeSet{i}.fastInner && d > timeInt.error(j)
                            % the distance d from the over-approximation to
                            % the halfspace is larger than the error, so the
                            % inner-approximation will definitely intersect
                            % the unsafe set -> provably unsafe in this
                            % special case
                            res = false;                     
                            if nargout == 2
                                R = aux_initReachableSet(RtimePoint,RtimeInt);
                            end
                            return
                        end
                        
                        % compute inner-approximation (unless already
                        % computed) and check intersection with unsafe set
                        [R_i,ind_Ri,ind_now] = aux_compInnerApprox(...
                            R_i,ind_Ri,timeInt.set,timeInt.error,j);
                        d = aux_intersectionCheck(R_i{ind_now},params.unsafeSet{i}.set);
                        
                        if d < 0
                            % intersection between inner-approximation and
                            % unsafe set -> provably unsafe
                            res = false; 
                            if nargout == 2
                                R = aux_initReachableSet(RtimePoint,RtimeInt);
                            end
                            return
                        elseif d < d_Fi
                            d_Fi = d; ind_Fi = j;  
                        end
                    else
                        % no intersection of current set timeInt.set{j}
                        % with i-th unsafe set -> extend already verified
                        % time intervals
                        params.unsafeSet{i}.time_ = ...
                            setdiff(params.unsafeSet{i}.time_,time(j,:));

                        % append to reachable set
                        allSAT(j,length(params.safeSet)+i) = true;
                        allSATbefore(j) = false;
                    end
                else
                    % append to reachable set
                    allSAT(j,length(params.safeSet)+i) = true;
                end
            end
        end

        % append sets which satisfy all specs to reachable set
        if nargout == 2
            idxSAT = all(allSAT,2) & ~allSATbefore;
            if any(idxSAT)
                [RtimePoint,RtimeInt] = aux_appendSets(RtimePoint,RtimeInt,...
                    timeInt.set,timePoint.set,timeInt.error,timePoint.error,time,idxSAT);
            end
        end
        
        % termination condition for do-until-loop
        if ((d_Go <= 0) && (d_Fo <= 0)) || (d_Gi > 0) || (d_Fi == 0)
        	break
        end
        
        % update max. allowed error
        if -d_Gi < d_Fi
            d = -d_Gi; ind = ind_Gi;
        else
            d = d_Fi; ind = ind_Fi;
        end
        
        if ~isempty(params.unsafeSet) && d_Fo ~= 0 && d_Fo < d
            d = d_Fo; ind = ind_Fo;
        end
        
        if ~isempty(params.safeSet) && d_Go >= 0 && d_Go < d
            d = d_Go; ind = ind_Go;
        end
        
        % update error value
        emax = max(0.1*emax,min(d*emax/timeInt.error(ind),0.9*emax));
    end
    
    % assign result
    res = (d_Go <= 0) && (d_Fo <= 0);

    % create reachSet object if desired
    if nargout == 2
        R = aux_initReachableSet(RtimePoint,RtimeInt);
    end
end


% Auxiliary functions -----------------------------------------------------

% safe/unsafe set manipulation
function [safeSet,unsafeSet] = aux_getSetsFromSpec(spec)
% extract polytopic safe sets and unsafe sets from the specifications; we
% also convert unsafe sets represented by halfspaces to safe sets and
% compute the interval enclosure for unsafe sets if they are bounded to
% speed up some intersection checks

    % all specifications must be of type 'safeSet' or 'unsafeSet'
    if ~all(arrayfun(@(s) any(strcmp(s.type,{'safeSet','unsafeSet'})),...
            spec,'UniformOutput',true))
        throw(CORAerror('CORA:notSupported',...
            'Only specification types ''safeSet'' and ''unsafeSet'' are supported.'));
    end

    % pre-allocate length, truncate at the end
    numSpec = numel(spec);
    safeSet = cell(numSpec,1); unsafeSet = cell(numSpec,1);
    idxSafeSet = 1; idxUnsafeSet = 1;
    
    for i=1:numel(spec)
        % convert set to a polytope and normalize constraint vectors
        P_norm = normalizeConstraints(polytope(spec(i).set),'A');
        
        if strcmp(spec(i).type,'safeSet')
            safeSet{idxSafeSet}.set = P_norm;
            safeSet{idxSafeSet}.time = spec(i).time;
            idxSafeSet = idxSafeSet + 1;

        elseif strcmp(spec(i).type,'unsafeSet')
            if representsa_(P_norm,'halfspace',eps)
                % if the unsafe set is only a single halfspace, convert to
                % a safe set using easy-to-compute set complement
                safeSet{idxSafeSet}.set = ~P_norm;
                safeSet{idxSafeSet}.time = spec(i).time;
                idxSafeSet = idxSafeSet + 1;
            else
                unsafeSet{idxUnsafeSet}.set = P_norm;
                unsafeSet{idxUnsafeSet}.time = spec(i).time;
                % pre-compute boundedness property: if it is bounded, we
                % use the interval enclosure for a quick intersection check
                unsafeSet{idxUnsafeSet}.isBounded = isBounded(P_norm);
                if unsafeSet{idxUnsafeSet}.isBounded
                    unsafeSet{idxUnsafeSet}.int = interval(P_norm);
                end
                idxUnsafeSet = idxUnsafeSet + 1;
            end
        end
    end

    % truncate cell arrays according to collected specs
    safeSet = safeSet(1:idxSafeSet-1);
    unsafeSet = unsafeSet(1:idxUnsafeSet-1);
end

% initializations
function err = aux_initError(sys,params)
% estimate the initial value for the maximum allowed error from simulations

    % initialize error value
    err = Inf;

    % initial points: center of facets from interval enclosure, where we
    % limit the number of points to a maximum of 20
    I = interval(params.R0); r = rad(I); Z = zonotope(I);
    
    if isempty(Z.G)
        points = Z.c;
    else
        if size(Z.G,2) <= 10
            points = [Z.c + Z.G, Z.c - Z.G];
        else
            [~,ind] = sort(r(r~=0),'descend');
            points = [Z.c + Z.G(:,ind(1:10)), Z.c - Z.G(:,ind(1:10))];
        end

        if ~representsa_(Z,'interval',eps)
            points = aux_getInitPoints(params.R0,points);
        end
    end
    
    % initialization
    paramSim.tFinal = params.tFinal;
    if isfield(params,'U')
        if isfield(params,'u')
            paramSim.u = params.u + center(params.U);
        else
            paramSim.u = center(params.U);
        end
    end
    
    % disturbance (required for getfcn... not yet added by validateOptions)
    paramSim.w = center(params.W);
    
    % loop over all initial points
    for i = 1:size(points,2)
        
        % simulate trajectories
        paramSim.x0 = points(:,i);
        if i == 1
            [t,~,~,y] = simulate(sys,paramSim);
            dt = paramSim.tFinal/ceil(paramSim.tFinal/(2*mean(diff(t))));
            sysDT = linearSysDT(sys,dt);
        else
            [t,~,~,y] = simulate(sysDT,paramSim);
        end
        
        % compute distance of output trajectory points to unsafe sets
        for k = 1:numel(params.unsafeSet)
            % read out k-th unsafe set
            P = params.unsafeSet{k}.set;

            for j = 1:size(y,1)
                % check whether unsafe set is active
                if contains(params.unsafeSet{k}.time,t(j))
                    scale = 10 - 9/t(end)*t;

                    % update error
                    if contains(P,y(j,:)')
                        % trajectory is contained in unsafe set
                        err_estimate = min(abs(P.A*y(j,:)'-P.b));
                        err = min(err,err_estimate*min(scale));

                    elseif params.unsafeSet{k}.isBounded
                        % trajectory outside of bounded unsafe set
                        % -> compute distance to interval enclosure
                        if aux_distanceIntervalPoint(params.unsafeSet{k}.int,y(j,:)') < err
                            err_estimate = aux_distancePolyPoint(P,y(j,:)');
                            err = min(err,err_estimate*min(scale));
                        end
                    else
                        % trajectory outside of unbounded unsafe set
                        % -> compute distance to unsafe set
                        if size(P.A,1) == 1
                            err_estimate = P.A*y(j,:)'-P.b;
                        else
                            err_estimate = aux_distancePolyPoint(P,y(j,:)');
                        end
                        err = min(err,err_estimate*min(scale));
                    end
                end
            end
        end
           
        % compute distance of output trajectory points to safe sets
        for k = 1:numel(params.safeSet)
            % read out k-th safe set
            P = params.safeSet{k}.set;

            if size(P.A,1) == 1 && representsa(params.safeSet{k}.time,"emptySet")
                scale = 10 - 9/t(end)*t;
                err_ = min(min(abs(P.A*y'-P.b),[],1) .* scale');
                err = min(err,err_);
            else
                % general case
                for j = 1:size(y,1)
                   if contains(params.safeSet{k}.time,t(j))
                       if contains(P,y(j,:)')
                           err_estimate = min(abs(P.A*y(j,:)'-P.b));
                       elseif size(P.A,1) == 1
                           err_estimate = P.A*y(j,:)'-P.b;
                       else
                           err_estimate = aux_distancePolyPoint(P,y(j,:)');
                       end
                       scale = 10 - 9/t(end)*t;
                       err = min(err,err_estimate*scale); 
                   end
                end
            end
        end        
    end
end

function points = aux_getInitPoints(R0,points)
% find 2n points on the boundary of the zonotope that are closest to the
% given points on the boundary of the interval enclosure

    % get zonotope properties
    c = R0.c; G = R0.G; [n,numGens] = size(G);

    % loop over all points
    for i = 1:size(points,2)
        
        % set-up linear program to minimize the norm 1 distance: 
        % min |points(:,i) - c + G*\alpha| s.t. -1 <= \alpha <= 1
        problem.f = [zeros(numGens,1); ones(2*n,1)];
        problem.Aeq = [-G, eye(n), -eye(n)];
        problem.beq = points(:,i) - c;
        problem.Aineq = [eye(numGens) zeros(numGens,2*n); ...
                        -eye(numGens) zeros(numGens,2*n); ...
                        zeros(2*n,numGens), -eye(2*n)];
        problem.bineq = [ones(2*numGens,1); zeros(2*n,1)];
        problem.lb = [];
        problem.ub = [];
        
        % solve linear program
        x = CORAlinprog(problem);
        points(:,i) = c + G*x(1:numGens);
    end
end

% inner-approximation
function [R_i,ind_Ri,ind_now] = aux_compInnerApprox(R_i,ind_Ri,R_o,Rout_error,ind_Ro)
% on-demand computation of inner-approximation based on outer-approximation
% save corresponding indices of Rout in ind_Ri
% return index corresponding to ind_Ro

% check if inner-approximation has already been computed:
% search for ind_Ro in ind_Ri list
ind_now = find(ind_Ri == ind_Ro,1,'first');

if ~isempty(ind_now)
    % inner-approximation is already computed, return same list + index
else
    % inner-approximation needs to be computed
    q = dim(R_o{ind_Ro});
    B = polytope((sqrt(q)*Rout_error(ind_Ro)*[eye(q),-eye(q)]')');
    R_i{end+1,1} = minkDiff(conZonotope(R_o{ind_Ro}),B);
    ind_Ri = [ind_Ri; ind_Ro];
    ind_now = length(ind_Ri);
end

end

% containment checks
function dist = aux_containmentCheckZono(Z,P)
% check if the zonotope Z is contained in the polytope P and return the
% distance

    % parameter for zonotope and polytope
    c = center(Z); G = generators(Z);
    C = P.A; d = P.b;
    
    % distance
    dist = max(-d + C*c + sum(abs(C*G),2)); 
end

function dist = aux_containmentCheckConZono(cZ,P)
% check if the constrained zonotope cZ is contained in the polytope P and
% return the distance

    % parameters for constrained zonotope and polytope
    c = cZ.c; G = cZ.G; A = cZ.A; b = cZ.b; l = size(G,2);
    C = P.A; d = P.b;
    
    % distance
    dist = -Inf;

    % init linprog struct
    problem.lb = -ones(l,1);
    problem.ub = ones(l,1);
    problem.Aeq = A;
    problem.beq = b;
    problem.Aineq = [];
    problem.bineq = [];
    
    for i = 1:length(d)
        problem.f = -C(i,:)*G;
        [~,dist_] = CORAlinprog(problem);
        dist = max(dist,C(i,:)*c - d(i) - dist_);
    end
end

% intersection checks
function dist = aux_intersectionCheck(cZ,P)
% check if the constrained zonotope cZ is intersecting the polyotpe P and 
% return the distance    

    cZ = conZonotope(cZ);

    % parameter for constrained zonotope and polytope
    c = cZ.c; G = cZ.G; A = cZ.A; b = cZ.b; l = size(G,2);
    C = P.A; d = P.b; n = size(C,2);

    % compute the distance by linear programming using the optimization
    % variables y = [a \alpha x] with slack variable a
    % 
    %   min 1^T a
    %
    %   s.t. C*x <= d, A*\alpha = b, -1 <= \alpha <= 1, 
    %        c + G*\alpha - x <= a, c + G*\alpha - x >= -a
    
    problem.f = [ones(n,1);zeros(n+l,1)];
    
    problem.Aineq = [zeros(size(C,1),n+l), C; ...
                     -eye(n), G, -eye(n); ...
                     -eye(n), -G, eye(n); ...
                     zeros(l,n), eye(l), zeros(l,n); ...
                     zeros(l,n), -eye(l), zeros(l,n)];
    problem.bineq = [d; -c; c; ones(l,1); ones(l,1)];
    
    problem.Aeq = [zeros(size(A,1),n) A zeros(size(A,1),n)];
    problem.beq = b;

    problem.lb = [];
    problem.ub = [];

    [~,dist] = CORAlinprog(problem);
    
end

function res = aux_intersectionCheckFast(Z,F)
% fast hierarchical check if a zonotope intersects a polytope

    res = false;

    % fast check using halfspace containment
    c = Z.c; G = Z.G; C = F.set.A; d = F.set.b;

    if max(C*c - sum(abs(C*G),2) - d) > 0
       return; 
    end
    
    % more accurate checks using interval enclosures and linear programming 
    if F.isBounded
       if isIntersecting(interval(Z),F.int) && ...
            	aux_intersectionCheck(Z,F.set) == 0 
            res = true;
       end
    elseif aux_intersectionCheck(Z,F.set) == 0
        res = true;
    end
end

% distance computations
function dist = aux_distanceIntersectionZono(Z,P)
% return maximum distance between a zonotope and polytope that intersect
    w = rad(interval(conZonotope(Z) & P));
    dist = 2*sqrt(sum(w.^2));
end

function d = aux_distanceIntervalPoint(I,p)
% compute the norm-1 distance between an interval and a point

    d = 0;
    l = infimum(I); u = supremum(I);
    
    for i = 1:length(p)
        if p(i) > u(i)
            d = d + u(i) - p(i); 
        elseif p(i) < l(i)
            d = d + l(i) - p(i);
        end
    end
end

function d = aux_distancePolyPoint(P,p)
% compute the norm-1 distance between a polytope P: C*x <= d and a point p

    % get polytope properties
    C = P.A; d = P.b;
    n = length(p); nrCon = size(C,1);
    
    % check how many halfspace constraints are violated
    offset = C*p - d;
    ind = find(offset > 0);
    
    if length(ind) == 1
        
        % only one halfspace constraint violated -> distance to polytope is
        % equal to the distance to the halfspace constraint
        d = offset(ind(1));
        
    elseif length(ind) == n
        
        % compute the vertex that is closest to the point by combining 
        % the violated halfspace constraints
        v = linsolve(C(ind,:),d(ind));
        d = sqrt(sum((v-p).^2));
        
    else
        % set-up linear program to minimize the norm 1 distance: 
        % min ||p - x||_1 s.t. C*x <= d
        problem.f = [zeros(n,1); ones(2*n,1)];
        problem.Aeq = [eye(n) eye(n) -eye(n)];
        problem.beq = p;
        problem.Aineq = [C zeros(nrCon,2*n); zeros(2*n,n) -eye(2*n)];
        problem.bineq = [d; zeros(2*n,1)];
        problem.lb = [];
        problem.ub = [];

        % solve linear program
        [~,d] = CORAlinprog(problem);
    end
end

function [RtimePoint,RtimeInt] = aux_appendSets(RtimePoint,RtimeInt,...
    R_o,R_o_tp,R_o_error,R_o_tp_error,time,idxSAT)
    
    % time-point solution
    % TODO: remove duplicates
    RtimePoint.set = [RtimePoint.set; R_o_tp(idxSAT)];
    RtimePoint.time = [RtimePoint.time; num2cell(time(idxSAT,2))];
    RtimePoint.error = [RtimePoint.error; R_o_tp_error([false;idxSAT])];

    % time-interval solution
    RtimeInt.set = [RtimeInt.set; R_o(idxSAT)];
    RtimeInt.error = [RtimeInt.error; R_o_error(idxSAT)];

    % cell-array for time-interval solution more cumbersome...
    RtimeInt_time = {};
    for jj=1:length(idxSAT)
        if idxSAT(jj)
            RtimeInt_time = [RtimeInt_time; {interval(time(jj,1),time(jj,2))}];
        end
    end
    RtimeInt.time = [RtimeInt.time; RtimeInt_time];

end

function [RtimePoint,RtimeInt] = aux_removeDuplicates(RtimePoint,RtimeInt)

% time-point solutions: if same time point, choose set with smaller error
time = cell2mat(RtimePoint.time);
i = 1;
while i < length(RtimePoint.error)
    % which indices have the same value as the i-th entry?
    idx = withinTol(time(i),time);
    % how many identical values?
    idx_no = nnz(idx);
    
    % go to next index if no identical values
    if idx_no > 1

        % look which entry has the smallest error value
        [~,minIdx] = min(RtimePoint.error(idx));
        idxNonzero = [1:length(idx)]' .* idx;
        idxNonzero = idxNonzero(idxNonzero > 0);
        idxRemoved = idxNonzero([1:minIdx-1,minIdx+1:end]);
        % remove entries from set, time, and errors
        RtimePoint.set(idxRemoved) = [];
        RtimePoint.time(idxRemoved) = [];
        RtimePoint.error(idxRemoved) = [];

        % update
        time = cell2mat(RtimePoint.time);

    end

    % increment counter
    i = i+1;

end

% time-interval solutions
i = 1;
time = zeros(length(RtimeInt.time),2);
for j=1:length(RtimeInt.time)
    time(j,:) = [infimum(RtimeInt.time{j}), supremum(RtimeInt.time{j})];
end

while i < length(RtimeInt.error)
    
    if any(all(time(i,1) >= time(i+1:end,1) & time(i,2) <= time(i+1:end,2),2))
        % remove as interval contained in larger time interval
        RtimeInt.set(i) = [];
        RtimeInt.time(i) = [];
        RtimeInt.error(i) = [];
        time(i,:) = [];
    else
        % increment counter
        i = i+1;
    end
end

end

function R = aux_initReachableSet(RtimePoint,RtimeInt)
    % remove duplicates (overlapping interval, multiple time points)
    [RtimePoint,RtimeInt] = aux_removeDuplicates(RtimePoint,RtimeInt);
    % splice reachSet object from different iterations in main loop
    R = reachSet(RtimePoint,RtimeInt);
    % re-order reachable set (note: time-interval solutions overlap)
    R = order(R);
end

% ------------------------------ END OF CODE ------------------------------
