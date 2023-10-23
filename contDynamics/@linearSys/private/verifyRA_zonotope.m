function [res,R] = verifyRA_zonotope(sys,params,options,spec)
% verifyRA_zonotope - verification for linear systems
%    via reach-avoid with zonotopes
%
% Syntax:
%    [res,R] = verifyRA_zonotope(sys,params,spec)
%
% Inputs:
%    sys - linearSys object
%    params - model parameters
%    spec - object of class specification
%
% Outputs:
%    res - true/false (true if specifications verified, otherwise false)
%    R - reachSet object with outer-approx. of the reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/verify

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       16-August-2021
% Last update:   08-November-2022 (MW, proper integration in master)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % get safe sets G and unsafe sets F
    [params.safeSet,params.unsafeSet] = aux_getSetsFromSpec(spec);
    % idea: allow to overwrite initial error (skipping initError call)

    % input argument pre-processing
    options = validateOptions(sys,mfilename,params,options);
    
    % estimate the required max. allowed error with simulations
    emax = aux_initError(sys,options,options.safeSet,options.unsafeSet);

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
            disp("Iteration " + cnt + ", current error: " + options.error);
        end
        [timeInt,timePoint,res,savedata] = reach_adaptive(sys,options);
        
        % save first time-point solution (same in every iteration)
        if cnt == 1
            RtimePoint.set = [RtimePoint.set; timePoint.set(1)];
            RtimePoint.time = [RtimePoint.time; {options.tStart}];
            RtimePoint.error = [RtimePoint.error; timePoint.error(1)];
        end
        
        % already falsified or (completely!) verified in reach_adaptive
        if ~res || ...
                ( (isempty(savedata.safeSet_unsat) || ...
                    isempty(cell2mat(savedata.safeSet_unsat))) && ...
                (isempty(savedata.unsafeSet_unsat) || ...
                    isempty(cell2mat(savedata.unsafeSet_unsat))) )
            % TODO: add missing reachable sets!
            return
        end
        
        % generate sequence of time intervals
        time = repmat(cell2mat(timePoint.time(2:end)),1,2);
        time(:,1) = [0; time(1:end-1,1)];
        % flag to see which sets can be saved in R
        allSAT = false(length(timeInt.set),length(options.safeSet)+length(options.unsafeSet));
        allSATbefore = true(length(timeInt.set),1);
        % save data for next run
        options.savedata = savedata;
        
        
        % clear variables for inner-approximation between runs
        R_i = {}; ind_Ri = [];
        
        % check if the reachable set is inside the safe sets       
        for i = 1:length(options.safeSet)
            for j = 1:length(timeInt.set)
                
                % check if specification is active and whether it has been
                % verified before
                if ( representsa(options.safeSet{i}.time,"emptySet") || ...
                    aux_intersectionCheckInterval(options.safeSet{i}.time,time(j,:)) ) ...
                        && aux_checkSet(options.savedata.safeSet_unsat{i},time(j,:))
                
                    % outer-approximation
                    d = aux_containmentCheckZono(timeInt.set{j},options.safeSet{i}.set);
                    
                    if d < 0
                        % outer-approximation is fully contained in the
                        % safe set -> extend already verified time intervals
                        options.savedata.safeSet_unsat{i} = ...
                            aux_removeFromUnsat(options.savedata.safeSet_unsat{i},time(j,:));

                        % append to reachable set
                        allSAT(j,i) = true;
                        allSATbefore(j) = false;
                    end
                    
                    if d > d_Go
                        d_Go = d; ind_Go = j; 
                    end
                    
                    % quick check for inner-approximation
                    if options.safeSet{i}.fastInner && d > timeInt.error(j)
                        % the distance d from the over-approximation to
                        % the safe set is larger than the contained
                        % error so that the inner-approximation will
                        % never be fully contained in the safe set
                        % -> provably unsafe in this special case
                        res = false; return;
                    end
                    
                    % inner-approximation
                    if d > timeInt.error(j)
                        % outer-approximation not fully contained in the
                        % safe set -> check whether inner-approximation
                        % intersects the safe set
                        [R_i,ind_Ri,ind_now] = aux_compInnerApprox(...
                            R_i,ind_Ri,timeInt.set,timeInt.error,j);
                        d = aux_containmentCheckConZono(R_i{ind_now},options.safeSet{i}.set);
                        if d > 0
                            % inner-approximation not fully contained in
                            % the safe set -> provably unsafe
                            res = false; return
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
        for i = 1:length(options.unsafeSet)
            for j = 1:length(timeInt.set)
                
                % check if specification is active (empty .time means that
                % the specification covers the entire time horizon)
                if ( representsa_(options.unsafeSet{i}.time,'emptySet') || ...
                        aux_intersectionCheckInterval(options.unsafeSet{i}.time,time(j,:)) ) ...
                        && aux_checkSet(options.savedata.unsafeSet_unsat{i},time(j,:))
                
                    % fast intersecting check for speed-up
                    if aux_intersectionCheckFast(timeInt.set{j},options.unsafeSet{i})

                        % outer-approximation
                        d = aux_distanceIntersectionZono(timeInt.set{j},options.unsafeSet{i}.set);

                        if d > d_Fo
                            d_Fo = d; ind_Fo = j; 
                        end

                        if options.unsafeSet{i}.fastInner && d > timeInt.error(j)
                            % the distance d from the over-approximation to
                            % the halfspace is larger than the error, so the
                            % inner-approximation will definitely intersect
                            % the unsafe set -> provably unsafe in this
                            % special case
                            res = false; return
                        end
                        
                        % compute inner-approximation (unless already
                        % computed) and check intersection with unsafe set
                        [R_i,ind_Ri,ind_now] = aux_compInnerApprox(...
                            R_i,ind_Ri,timeInt.set,timeInt.error,j);
                        d = aux_intersectionCheck(R_i{ind_now},options.unsafeSet{i}.set);
                        
                        if d < 0
                            % intersection between inner-approximation and
                            % unsafe set -> provably unsafe
                            res = false; return
                        elseif d < d_Fi
                            d_Fi = d; ind_Fi = j;  
                        end
                    else
                        % no intersection of current set timeInt.set{j} with F{i}
                        % -> extend already verified time intervals
                        options.savedata.unsafeSet_unsat{i} = ...
                            aux_removeFromUnsat(options.savedata.unsafeSet_unsat{i},time(j,:));

                        % append to reachable set
                        allSAT(j,length(options.safeSet)+i) = true;
                        allSATbefore(j) = false;
                    end
                else
                    % append to reachable set
                    allSAT(j,length(options.safeSet)+i) = true;
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
        
        if ~isempty(options.unsafeSet) && d_Fo ~= 0 && d_Fo < d
            d = d_Fo; ind = ind_Fo;
        end
        
        if ~isempty(options.safeSet) && d_Go >= 0 && d_Go < d
            d = d_Go; ind = ind_Go;
        end
        
        % update error value
        emax = max(0.1*emax,min(d*emax/timeInt.error(ind),0.9*emax));
    end
    
    % assign result
    res = (d_Go <= 0) && (d_Fo <= 0);

    % create reachSet object if desired
    if nargout == 2
        % remove duplicates (overlapping interval, multiple time points)
        [RtimePoint,RtimeInt] = aux_removeDuplicates(RtimePoint,RtimeInt);
        % splice reachSet object from different iterations in main loop
        R = reachSet(RtimePoint,RtimeInt);
        % re-order reachable set (note: time-interval solutions overlap)
        R = order(R);
    end
end


% Auxiliary functions -----------------------------------------------------

% safe/unsafe set manipulation
function [G,F] = aux_getSetsFromSpec(spec)
% extract safe sets G and unsafe sets F from the specifications

    G = {}; F = {};
    
    for i = 1:length(spec)
        if strcmp(spec(i).type,'safeSet')
            G{end+1}.set = normalizeConstraints(polytope(spec(i).set),'A');
            G{end}.time = spec(i).time;
        elseif strcmp(spec(i).type,'unsafeSet')
            tmp = normalizeConstraints(polytope(spec(i).set),'A');
            if size(tmp.A,1) > 1
                F{end+1}.set = tmp;
                F{end}.time = spec(i).time;
                if size(F{end}.set.A,1) > size(F{end}.set.A,2)
                   F{end}.int = interval(F{end}.set);
                   F{end}.isBounded = ~(any(isinf(infimum(F{end}.int))) ...
                                       | any(isinf(supremum(F{end}.int))));
                else
                   F{end}.isBounded = false; 
                end
            else
                G{end+1}.set = polytope(-tmp.A,-tmp.b);
                G{end}.time = spec(i).time;
            end
        else
            throw(CORAerror('CORA:notSupported',...
                'Given type of specification is not supported.'));
        end
    end
end

% initializations
function err = aux_initError(sys,options,G,F)
% estimate the initial value for the maximum allowed error from simulations

    % initial points: center of facets from interval enclosure, where we
    % limit the number of points to a maximum of 20
    I = interval(options.R0); r = rad(I);
    Z = zonotope(I); c = center(Z); g = generators(Z);
    
    if ~isempty(g)
        if size(g,2) <= 10
            points = [c + g, c - g];
        else
            [~,ind] = sort(r(r~=0),'descend');
            points = [c + g(:,ind(1:10)), c - g(:,ind(1:10))];
        end

        if ~representsa_(Z,'interval',eps)
            points = aux_getInitPoints(options.R0,points);
        end
    else
        points = c;
    end
    
    % initialization
    paramSim.tFinal = options.tFinal;
    if isfield(options,'U')
        if isfield(options,'u')
            paramSim.u = options.u + center(options.U);
        else
            paramSim.u = center(options.U);
        end
    end
    
    % disturbance (required for getfcn... not yet added by validateOptions)
    paramSim.w = center(options.W);
        
    err = Inf;
    
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
        for k = 1:length(F)
            for j = 1:size(y,1)
               if contains(F{k}.time,t(j))
                   scale = 10 - 9/t(end)*t;
                   if contains(F{k}.set,y(j,:)')
                       temp = min(abs(F{k}.set.A*y(j,:)'-F{k}.set.b));
                       err = min(err,temp*scale);
                   elseif F{k}.isBounded 
                       if aux_distanceIntervalPoint(F{k}.int,y(j,:)') < err
                           temp = aux_distancePolyPoint(F{k}.set,y(j,:)');
                           err = min(err,temp*scale);
                       end
                   else
                       if size(F{k}.set.A,1) == 1
                           temp = F{k}.set.A*y(j,:)'-F{k}.set.b;
                       else
                           temp = aux_distancePolyPoint(F{k}.set,y(j,:)');
                       end
                       err = min(err,temp*scale);
                   end
               end
            end
        end
           
        % compute distance of output trajectory points to safe sets
        for k = 1:length(G)
            if size(G{k}.set.A,1) == 1 && representsa(G{k}.time,"emptySet")
                scale = 10 - 9/t(end)*t;
                tmp = min(abs(G{k}.set.A*y'-G{k}.set.b),[],1);
                err_ = min(tmp .* scale');
                err= min(err,err_);
            else
                for j = 1:size(y,1)
                   if contains(G{k}.time,t(j))
                       if contains(G{k}.set,y(j,:)')
                           temp = min(abs(G{k}.set.A*y(j,:)'-G{k}.set.b));
                       elseif size(G{k}.set.A,1) == 1
                           temp = G{k}.set.A*y(j,:)'-G{k}.set.b;
                       else
                           temp = aux_distancePolyPoint(G{k}.set,y(j,:)');
                       end
                       scale = 10 - 9/t(end)*t;
                       err = min(err,temp*scale); 
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
    c = center(R0); G = generators(R0); [n,m] = size(G);

    % loop over all points
    for i = 1:size(points,2)
        
        % set-up linear program to minimize the norm 1 distance: 
        % min |points(:,i) - c + G*\alpha| s.t. -1 <= \alpha <= 1
        f = [zeros(m,1); ones(2*n,1)];
        Aeq = [-G eye(n) -eye(n)]; beq = points(:,i) - c;
        A = [eye(m) zeros(m,2*n); -eye(m) zeros(m,2*n); ...
             zeros(2*n,m), -eye(2*n)];
        b = [ones(2*m,1); zeros(2*n,1)];
        
        % solve linear program
        options = optimoptions('linprog','Display','off');
        x = linprog(f,A,b,Aeq,beq,[],[],options);
        points(:,i) = c + G*x(1:m);
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
    B = polytope(sqrt(q)*Rout_error(ind_Ro)* [eye(q),-eye(q)]');
    R_i{end+1,1} = minkDiff(conZonotope(R_o{ind_Ro}),B);
    ind_Ri = [ind_Ri; ind_Ro];
    ind_now = length(ind_Ri);
end

% previously: full computation of inner-approximation of the reachable set
% R_i = cell(size(R_o));
% for i = 1:length(R_o)
%    B = polytope(sqrt(n)*actError.timeInt(i)* [eye(n),-eye(n)]');
%    R_i{i} = conZonotope(R_o{i}) - B;
% end

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

    % parameter for constrained zonotope and polytope
    c = cZ.c; G = cZ.G; A = cZ.A; b = cZ.b; l = size(G,2);
    C = P.A; d = P.b;
    
    % distance
    options = optimoptions('linprog','display','off');
    dist = -inf;
    
    for i = 1:length(d)
        [~,dist_] = linprog(-C(i,:)*G,[],[],A,b,-ones(l,1),ones(l,1), ...
                                                           [],options);                   
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
    
    f = [ones(n,1);zeros(n+l,1)];
    
    Ain = [zeros(size(C,1),n+l),C; -eye(n) G -eye(n); -eye(n) -G eye(n); ...
           zeros(l,n) eye(l) zeros(l,n); zeros(l,n) -eye(l) zeros(l,n)];
    bin = [d; -c; c; ones(l,1); ones(l,1)];
    
    Aeq = [zeros(size(A,1),n) A zeros(size(A,1),n)];
    beq = b;
    
    linprogoptions = optimoptions('linprog','display','off');
    
    [~,dist] = linprog(f,Ain,bin,Aeq,beq,[],[],[],linprogoptions);
end

function res = aux_intersectionCheckFast(Z,F)
% fast hierarchical check if a zonotope intersects a polytope

    res = false;

    % fast check using halfspace containment
    c = center(Z); G = generators(Z);
    C = F.set.A; d = F.set.b;

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

function res = aux_intersectionCheckInterval(I1,I2)
% copied and shortened from interval/isIntersecting/isIntersecting1D to
% avoid costly instantiation of interval objects (I1 is still an interval)

    % read infimum and supremum
    I1 = [infimum(I1), supremum(I1)];

    res = false;
    if I1(1) <= I2(1)
        if I2(1) <= I1(2)
            res = true;
        end
    else % I2(1) < I1(1)
        if I(1) <= I2(2)
            res = true;
        end
    end

end

% distance computations
function dist = aux_distanceIntersectionZono(Z,P)
% return maximum distance between a zonotope and polytope in the case that
% both are intersecting

    w = rad(interval(conZonotope(Z) & P));
    dist = 2*sqrt(sum(w.^2));
end

function d = aux_distanceIntervalPoint(I,p)
% compute the norm 1 distance between a point p and an interval int

    d = 0;
    l = infimum(I);
    u = supremum(I);
    
    for i = 1:length(p)
       if p(i) > u(i)
          d = d + u(i) - p(i); 
       elseif p(i) < l(i)
          d = d + l(i) - p(i);
       end
    end
end

function d = aux_distancePolyPoint(P,p)
% compute the norm 1 distance between a point p and a polytope P: C*x <= d    

    % get polytope properties
    C = P.A; d = P.b;
    n = length(p); m = size(C,1);
    
    % check how many halfspace constraints are violated
    temp = C*p - d;
    ind = find(temp > 0);
    
    if length(ind) == 1
        
        % only one halfspace constraint violated -> distance to polytope is
        % equal to the distance to the halfspace constraint
        d = temp(ind(1));
        
    elseif length(ind) == n
        
        % compute the vertex that is closest to the point by combining 
        % the violated halfspace constraints
        v = linsolve(C(ind,:),d(ind));
        d = sqrt(sum((v-p).^2));
        
    else
        % set-up linear program to minimize the norm 1 distance: 
        % min ||p - x||_1 s.t. C*x <= d
        f = [zeros(n,1); ones(2*n,1)];
        Aeq = [eye(n) eye(n) -eye(n)]; beq = p;
        A = [C zeros(m,2*n); zeros(2*n,n) -eye(2*n)];
        b = [d; zeros(2*n,1)];

        % solve linear program
        options = optimoptions('linprog','Display','off');
        [~,d] = linprog(f,A,b,Aeq,beq,[],[],options);
    end
end

% safe/unsafe set time interval check
function doCheck = aux_checkSet(set_unsat,t)
% returns whether the current safe/unsafe set has to be checked for the 
% given time interval, which is not required if the verification has been
% successful in a prior iteration

doCheck = true;

% quick check if set is already fully verified
if isempty(set_unsat)
    doCheck = false;
    return;
end

% check for intersection
if t(2) <= set_unsat(1,1) || t(1) > set_unsat(end,2)
    % upper bound smaller than lower bound of first time interval
    % lower bound larger than upper bound of last time interval
    doCheck = false;
    return;
elseif any(t(1) >= set_unsat(1:end-1,2) & t(2) <= set_unsat(2:end,1))
    % bounds between unverified time intervals
    doCheck = false;
    return;
end

end

function set_unsat = aux_removeFromUnsat(set_unsat,t)
% adapt FGunsat so that timeInterval is not part of time intervals covered

FGunsat_col = reshape(set_unsat',numel(set_unsat),1);
if mod(sum(t(1) >= FGunsat_col),2) ~= 0
    % lower bound starts inside unverified time interval
    % t(1) \in [ timeInterval )
    idx = find(t(1) >= set_unsat(:,1) & t(1) <= set_unsat(:,2));
    
    if t(2) <= set_unsat(idx,2)
        if t(1) > set_unsat(idx,1)
            set_unsat = [set_unsat(1:idx-1,:); [set_unsat(idx,1), t(1)]; set_unsat(idx:end,:)];
            idx = idx + 1;
        end
        % split, potential merge later
        set_unsat(idx,1) = t(2);
        t = [];
    else
        % remove interval, potential merge later
        set_unsat(idx,2) = t(1);
        if idx < size(set_unsat,1)
            t(1) = set_unsat(idx+1,1);
        end
        if t(2) <= t(1)
            t = [];
        end
    end
end

% now: lower bound starts in between unverified time intervals or at
% maximum at the start point of an unverified set
% t(1) \in [ notTimeInterval )
while ~isempty(t)
    idx = find(t(1) <= set_unsat(:,1),1,'first');
    % upper bound is at least at the beginning of the next time interval
    if t(2) <= set_unsat(idx,2)
        % split, potential merge later
        set_unsat(idx,1) = t(2);
        t = [];
    else
        % remove entire thing (full time interval verified)
        if idx < size(set_unsat,1)
            t(1) = set_unsat(idx,2);
            if t(2) < set_unsat(idx+1,1)
                t = [];
            end
        else
            t = [];
        end
        set_unsat(idx,:) = [];
    end
end

% remove
idxRemove = abs(set_unsat(:,2) - set_unsat(:,1)) < 1e-14;
set_unsat(idxRemove,:) = [];

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
    temp = {};
    for jj=1:length(idxSAT)
        if idxSAT(jj)
            temp = [temp; {interval(time(jj,1),time(jj,2))}];
        end
    end
    RtimeInt.time = [RtimeInt.time; temp];

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

% ------------------------------ END OF CODE ------------------------------
