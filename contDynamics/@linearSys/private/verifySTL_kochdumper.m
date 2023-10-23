function [res,R,fals] = verifySTL_kochdumper(sys,params,options,spec)
% verifySTL_kochdumper - verification of linear systems vs. temporal logic specs [1]
%
% Syntax:
%    [res,R,fals] = verifySTL_kochdumper(sys,params,spec)
%
% Inputs:
%    sys - linearSys object
%    params - model parameters
%    spec - object of class specification
%
% Outputs:
%    res - boolean (true if specifications verified, otherwise false)
%    R - outer-approx. of the reachable set (class reachSet)
%    fals - struct storing the initial state and inputs for the falsifying
%           trajectory
%
% References:
%    [1] N. Kochdumper et al. "Fully-Automated Verification of Linear 
%        Time-Invariant Systems against Signal Temporal Logic 
%        Specifications via Reachability Analysis"
%    [2] H. Roehm et al. "STL Model Checking of Continuous and Hybrid
%        Systems", International Symposium on Automated Technology for 
%        Verification and Analysis, pp. 412-427, 2016.
%    [3] M. Wetzlinger et al. "Fully-Automated Verification of Linear 
%        Systems Using Inner- and Outer-Approximations of Reachable Sets"
%    [4] M. Althoff. "Reachability Analysis and its Application to the
%        Safety Assessment of Autonomous Cars", Dissertation 2010
%    [5] N. Kochdumper et al. "Provably Safe Reinforcement Learning via 
%        Action Projection using Reachability Analysis and Polynomial 
%        Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/verify

% Authors:       Niklas Kochdumper
% Written:       17-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input argument pre-processing
options = validateOptions(sys,mfilename,params,options);

% determine the required time horizon
if ~isfield(options,'u')
    options.tFinal = maximumTime(spec.set);
end

% initialization
pred.point = []; pred.int = []; predNeg.point = []; predNeg.int = [];

% choose initial time step size
[dt,options] = aux_initTimeStepSize(sys,options);

% preprocess temporal logic specifications
[phi,sets] = aux_preprocessTemporalLogic(spec.set);
[phiNeg,setsNeg] = aux_preprocessTemporalLogic(~spec.set);

% loop until specification can be proven or disproven
while true

    % compute reachable set
    [R,ind] = aux_compReachSet(sys,options,dt);

    % check if specification is satisfied
    [res,~,elem,pred] = aux_modelCheckReachSet(phi,R,ind,dt,sets,pred,false);

    if res
        fals = []; return;
    end

    % check if specification is violated (by checking negated formula)
    if ~elem

        [res,alpha,~,predNeg] = aux_modelCheckReachSet(phiNeg,R,ind,dt, ...
                                                setsNeg,predNeg,true);
    
        if ~res
            fals = aux_falsifyingTrajectory(alpha,options.R0,options.U,dt);
            R = []; return;
        end
    end

    % update time step size
    dt = dt/2;

    if isfield(options,'u')
        options.u = repelem(options.u,1,2);
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [R,ind] = aux_compReachSet(sys,options,dt)
% compute an over-approximation of the reachable set

    order = 10;
    n = dim(options.R0);
    p_x = size(generators(options.R0),2);
    p_u = size(generators(options.U),2);
    t = 0:dt:(options.tFinal + dt);

    % propagation matrices
    eAdt = expm(sys.A*dt);
    T = aux_constInputPropMat(sys.A,dt);

    % modified input sets
    u = sys.B * center(options.U);
    U = sys.B * (options.U + (-center(options.U)));

    if ~isempty(sys.c)
        u = u + sys.c;
    end

    if ~isfield(options,'u')
        options.u = zeros(n,length(t)-1);
    else
        options.u = sys.B*[options.u,options.u(:,end)];
    end
    
    % auxiliary terms
    E = aux_remMatrixExp(sys.A,dt,order);
    [Usum,Asum] = aux_precomputeAuxiliaryTerms(sys.A,U,dt,order);
    [F,G] = aux_curvatureEnclosure(sys.A,E,dt,order);

    % initialization
    Rint.set = cell(length(t)-1,1); Rint.time = cell(length(t)-1,1);
    Rpoint.set = cell(length(t),1); Rpoint.time = cell(length(t),1);

    ind.point = cell(length(t),1); ind.int = cell(length(t)-1,1);
    ind.point{1} = 1:p_x;

    H = cell(length(t),1); H{1} = zonotope(options.R0);
    P_u = zeros(n,1); P_U = zonotope(P_u); 
    eAt = eye(n); err = 0;
    C = F*options.R0 + G*u;

    Rpoint.time{1} = 0; Rpoint.set{1} = H{1};

    % loop over all time steps
    for i = 1:length(t)-1

         eAt = eAdt * eAt;

        H{i+1} = eAdt * H{i} + T*(u + options.u(:,i));

        P_U = eAdt * P_U + T*U;

        err = err + aux_errorParticularSolution(eAt,U,Usum,Asum,E,dt);

        Rpoint.set{i+1,1} = H{i+1} + P_U + err;
        Rint.set{i,1} = enclose(H{i},H{i+1}) + P_U + (interval(C) + err); 

        ind.point{i+1,1} = 1:p_x+i*p_u;
        ind.int{i,1} = [1:p_x,2*p_x+2:2*p_x+i*p_u+1]; 
        
        Rpoint.time{i+1,1} = t(i+1);
        Rint.time{i,1} = interval(t(i),t(i+1));

        C = eAdt * C;
    end

    % construct reachSet object
    R = reachSet(Rpoint,Rint);
end

function [res,alpha,elem,pred] = aux_modelCheckReachSet(phi,R,ind,dt,sets,pred,neg)
% check if the reachable set satisfies a signal temporal logic formula

    % initialization
    alpha = []; elem = false; 
    predInt = 2*ones(length(sets),length(R.timeInterval.set));
    predPoint = 2*ones(length(sets),length(R.timePoint.set));

    % reuse already verified times from previous iteration
    if ~isempty(pred.int)
        predInt(:,1:end-1) = repelem(pred.int(:,1:end-1),1,2);
    end

    % evaluate predicates on the time interval reachable set
    pred.int = aux_evaluatePredicates(R.timeInterval.set,sets,ind.int,predInt);

    % reuse already verified times from previous iteration
    if ~isempty(pred.point)
        predPoint(:,1:2:size(predPoint,2)) = pred.point(:,1:end-1);
    end

    % time point dosn't need to be checked if time interval set is verified
    for i = 1:size(pred.int,1)
        for j = 1:size(pred.int,2)
            if pred.int(i,j) < 2
                predPoint(i,j) = pred.int(i,j);
                predPoint(i,j+1) = pred.int(i,j);
            end
        end
    end

    % evaluate predicates on the time point reachable set
    pred.point = aux_evaluatePredicates(R.timePoint.set,sets,ind.point,predPoint);

    predInt = pred.int; predPoint = pred.point;
    predInt(predInt == 3) = 0; predPoint(predPoint == 3) = 0;
    pred.int(pred.int == 3) = 2; pred.point(pred.point == 3) = 2;

    % convert to sampled-time STL according to Sec. 4.2 in [2]
    phi = sampledTime(phi,dt,true,predInt,predPoint);
    phi = combineNext(phi);

    % catch the case where the whole formula is already true/false
    if strcmp(phi.type,'true')
        if ~neg
            res = true; return;
        else
            res = false; alpha=-1+2*rand(length(ind.point{end}),1); return;
        end
    elseif strcmp(phi.type,'false')
        if ~neg
            res = false; elem = true; return;
        else
            res = true; return;
        end
    end

    % check if conversion to conjunctive normal form brings a benefit
    if aux_benefitCNF(predPoint,predInt)
        
        % convert to Reachset Temporal Logic according to Sec. 4.1 in [2]
        [~,list] = stl2rtl(phi,[],true);

        % model check the formula in conjunctive normal form
        cons = aux_modelCheckCNF(list,R,ind,sets,neg);

    else

        % model check the original formula
        con = aux_precomputeConstraints(phi,R,ind,sets);
        cons = aux_modelCheck(phi,con,0);
    end

    % assign verification result
    if isempty(cons)
        res = true; return;
    else
        res = false;
    end

    % try to find a falsifying trajectory if the negated formula is checked
    if neg
        [res,alpha] = aux_falsifyingFactors(cons);
    end
end
    
function cons = aux_modelCheckCNF(list,R,ind,sets,neg)

    % check if reachable set satisfies the temporal logic formula
    cons = {};

    for i = 1:length(list)              % loop over all conjunctions

        resTmp = false; cons_= {};

        for j = 1:length(list{i})       % loop over all disjunctions

            % get constraints on factors from intersection between the
            % reachable set and the unsafe sets defined by logic formula
            eq = list{i}{j}.lhs.lhs;
            time = list{i}{j}.time;

            conTemp = aux_getConstraints(eq,R,ind,time,sets);

            % terminate loop as soon as one disjunction is satisfied
            if isempty(conTemp)
                resTmp = true; break;
            else
                cons_{end+1} = conTemp;
            end
        end

        % terminate loop as soon as one conjunction is not satisfied
        if ~resTmp

            cons_ = aux_orConstraints(cons_);

            if ~isempty(cons_)
                if ~neg
                    cons = cons_; return;
                elseif iscell(cons_)
                    cons = [cons,cons_];
                else
                    cons = []; return;
                end
            end
        end
    end
end

function con = aux_modelCheck(phi,cons,time)
% recursive function to compute the constraints on the zonotope factors
% that have to hold such that the temporal logic formula is satisfied

    if strcmp(phi.type,'next')

        con = aux_modelCheck(phi.lhs,cons,time + phi.from);

    elseif strcmp(phi.type,'globally') || strcmp(phi.type,'finally')

        con = aux_modelCheck(phi.lhs,cons,time + 0.5);

    elseif isempty(phi.id) && strcmp(phi.type,'|')

        con1 = aux_modelCheck(phi.lhs,cons,time);
        con2 = aux_modelCheck(phi.rhs,cons,time);

        list = {};

        for i = 1:length(con1)
            for j = 1:length(con2)
                list{end+1} = aux_combineConstraints(con1{i},con2{j});
            end
        end

        con = aux_andConstraints(list);

        for i = 1:length(con)
            if aux_emptyConstraint(con{i}.C,con{i}.d)
                con{i} = [];
            end
        end
        
        con = con(~cellfun('isempty',con));

    elseif isempty(phi.id) && strcmp(phi.type,'&')

        con1 = aux_modelCheck(phi.lhs,cons,time);
        con2 = aux_modelCheck(phi.rhs,cons,time);

        con = aux_andConstraints([con1,con2]);

    elseif ~phi.temporal

        if mod(time,1) == 0
            con = cons.timePoint{phi.id}{time+1};
        else
            con = cons.timeInt{phi.id}{floor(time)+1};
        end
    end
end

function res = aux_benefitCNF(predPoint,predInt)
% check if the conversion to conjunctive normal form brings a benefit

    res = true;

    % only one predicate -> no benefit since no disjunction of predicates
    if size(predPoint,1) == 1
        res = false; return;
    end

    res = false;
end

function conTemp = aux_getConstraints(eq,R,ind,time,sets)
% compute constraints on factors from intersection betwee the reachable set 
% and unsafe sets defined by the logic predicates

    % represent logic equation as a union of unsafe sets
    if isempty(eq.id)

        % convert logic equation to union of safe sets
        eq = disjunctiveNormalForm(eq);
        safeSet = getClauses(eq,'dnf');

        for k = 1:length(safeSet)
            if isempty(safeSet{k}.id) || iscell(sets{safeSet{k}.id})
                safeSet{k} = convert2set(safeSet{k});
            else
                safeSet{k} = sets{safeSet{k}.id};
            end
        end

        % convert to a union of unsafe sets
        unsafeSet = aux_safe2unsafe(safeSet);

    else
        
        % use precomputed set
        unsafeSet = sets{eq.id};

        if ~iscell(unsafeSet)
            unsafeSet = aux_safe2unsafe(unsafeSet);
        end
    end

    % select the correct reachable set
    if mod(time,1) == 0         % time point reachable set
        set = R.timePoint.set{time+1};
        index = ind.point{time+1};
    else                        % time interval reachable set
        set = R.timeInterval.set{floor(time)+1};
        index = ind.int{floor(time)+1};
    end

    % compute constraints on factors
    conTemp = aux_getFactaux_orConstraints(set,unsafeSet,index);
end

function cons = aux_getFactaux_orConstraints(set,unsafeSet,index)
% get constraints on the factors from intersection with unsafe sets

    cons = {};

    % loop over all unsafe sets
    for k = 1:length(unsafeSet)

        % check if there is an intersection with the unsafe set
        if isIntersecting(unsafeSet{k},set)

            % construct constraints
            if isa(unsafeSet{k},'halfspace')

                [C,d] = aux_getConstraintHalfspace(set,unsafeSet{k},index);
                cons{end+1}.C = C; cons{end}.d = d;
                    
            elseif contains(unsafeSet{k},set)

                % create fake constraint containing whole domain
                C = zeros(1,length(index)); C(1) = 1; d = 2;
                cons{end+1}.C = C; cons{end}.d = d;

            else

                % loop over all polytope halfspaces
                A = unsafeSet{k}.A; b = unsafeSet{k}.b;
                C = []; d = [];

                for i = 1:size(A,1)

                    ind = setdiff(1:size(A,1),i);
                    facet = conHyperplane(A(i,:),b(i),A(ind,:),b(ind));

                    if isIntersecting(facet,set)
                        [C_,d_] = aux_getConstraintHalfspace(set,facet.h,index);
                        C = [C; C_]; d = [d; d_];
                    end
                end

                cons{end+1}.C = C; cons{end}.d = d;
            end
        end
    end
end

function [C,d] = aux_getConstraintHalfspace(set,unsafeSet,index)
% get the constraint on the factors for an intersection with a halfspace

    % project set onto the halfspace normal direction
    C = unsafeSet.c' * generators(set);

    % get indices of factors that are not initial states or inputs
    ind = setdiff(1:size(C,2),index);

    % construct constraint offset and vector
    d = unsafeSet.d - unsafeSet.c' * center(set) + sum(abs(C(ind)));
    C = C(index);
end

function con = aux_precomputeConstraints(phi,R,index,sets)
% precompute the constraints on the factors for all predicates

    cnt = 1;

    % convert sets to unsafe sets
    for i = 1:length(sets)
        if ~iscell(sets{i})
            sets{i} = aux_safe2unsafe(sets{i});
        end
    end

    % determine times where constraints are active
    [timeInt,timePoint] = getTimes(phi);

    % loop over all constraints
    for i = 1:length(timeInt)

        % get constraints for time interval reachable sets
        if ~isempty(timeInt{i})

            ind = unique(cellfun(@(x) infimum(x),timeInt{i})) + 1;
    
            con.timeInt{i} = cell(max(ind),1);
    
            for j = 1:length(ind)
                con.timeInt{i}{ind(j)} = aux_getFactaux_orConstraints( ...
                     R.timeInterval.set{ind(j)},sets{i},index.int{ind(j)});
    
                for k = 1:length(con.timeInt{i}{ind(j)})
                    con.timeInt{i}{ind(j)}{k}.id = cnt; cnt = cnt + 1;
                end
            end
        end

        % get constraints for time point reachable sets
        if ~isempty(timePoint{i})

            ind = unique(cell2mat(timePoint{i})) + 1;
    
            con.timePoint{i} = cell(max(ind),1);
    
            for j = 1:length(ind)
                con.timePoint{i}{ind(j)} = aux_getFactaux_orConstraints( ...
                      R.timePoint.set{ind(j)},sets{i},index.point{ind(j)});
    
                for k = 1:length(con.timePoint{i}{ind(j)})
                    con.timePoint{i}{ind(j)}{k}.id = cnt; cnt = cnt + 1;
                end
            end
        end
    end
end

function list = aux_orConstraints(cons)
% or-connection of all constraints in the list

    % filter out constraints that contain the whole factor domain
    for i = 1:length(cons)
        for k = 1:length(cons{i})
        
            n = size(cons{i}{k}.C,2);
            dom = interval(-ones(n,1),ones(n,1));
    
            if contains(polytope(cons{i}{k}.C,cons{i}{k}.d),dom)
                cons{i} = []; break;
            end
        end
    end

    cons = cons(~cellfun('isempty',cons));

    if isempty(cons)
        list = 0; return;
    end

    % combine all constraints in the list by computing their intersection
    list = cons{1};

    for i = 2:length(cons)

        list_ = {};

        for k = 1:length(cons{i})
            for j = 1:length(list)

                con = aux_combineConstraints(cons{i}{k},list{j});

                if ~aux_emptyConstraint(con.C,con.d)
                    list_{end+1} = con;
                end
            end
        end

        list = list_;
    end
end

function list = aux_andConstraints(list)
% and-connection of all constraints in the list

    for i = 1:length(list)
        redundant = false;
        for j = 1:length(list)
            if i ~= j && ~isempty(list{j})
                id = intersect(list{i}.id,list{j}.id);
                if length(id) == length(list{j}.id)
                    redundant = true; break;
                end
            end
        end
        if redundant
            list{i} = [];
        end
    end

    if ~isempty(list)
        list = list(~cellfun('isempty',list));
    end
end

function con = aux_combineConstraints(con1,con2)
% combine two inequality constraints C*x < d

    p1 = size(con1.C,2); p2 = size(con2.C,2);

    if p1 > p2
        con.C = [con1.C; [con2.C, zeros(size(con2.C,1), p1 - p2)]];
    else
        con.C = [[con1.C, zeros(size(con1.C,1), p2 - p1)]; con2.C];
    end

    con.d = [con1.d;con2.d]; 

    if isfield(con1,'id') && isfield(con2,'id')
        id = [con1.id;con2.id];
        [con.id,ind] = unique(id);
        con.C = con.C(ind,:); con.d = con.d(ind,:);
    end
end

function res = aux_emptyConstraint(C,d)
% check if a constraint C*a < d is empty on the domain a \in [-1,1]

    n = size(C,2);

    if size(C,1) == 1

        dom = interval(-ones(n,1),ones(n,1));
        res = ~isIntersecting(halfspace(C,d),dom);

    else
        
        ub = ones(n,1); lb = -ones(n,1); f = ones(n,1);
        options = optimoptions('linprog','display','off');
    
        [~,~,exitflag] = linprog(f,C,d,[],[],lb,ub,options);
    
        res = exitflag ~= 1;
    end
end

function fals = aux_falsifyingTrajectory(alpha,R0,U,dt)
% construct a falsifying trajectory out of the factor values

    % compute initial state
    c_x = center(R0); G_x = generators(R0);
    fals.x0 = c_x + G_x * alpha(1:size(G_x,2));

    % compute input trajectory
    c_u = center(U); G_u = generators(U); m = size(G_u,2);
    fals.u = []; cnt = size(G_x,2) + 1;

    while cnt < length(alpha)
        fals.u = [fals.u c_u + G_u * alpha(cnt:cnt+m-1)];
        cnt = cnt + m;
    end

    % compute times
    fals.t = 0:dt:size(fals.u,2)*dt;
end

function [res,alpha] = aux_falsifyingFactors(cons)
% try to find factor values that falsify the temporal logic equation

    % determine number of variables
    n = 0;

    for i = 1:length(cons)
        n = max(n,size(cons{i}.C,2));
    end

    % initialization
    A = []; b = []; Aeq1 = []; Aeq2 = []; beq = []; intcon = [];
    B = polytope(interval(-ones(n,1),ones(n,1))); 

    for i = 1:length(cons)

        Aeq_ = []; C = cons{i}.C; d = cons{i}.d;

        for j = 1:size(C,1)

            % adapt dimension of constraint matrices
            if size(C,2) < n
                C = [C, zeros(size(C,1),n-size(C,2))];
            end

            % 1st constraint in Equation (9) in [5]
            temp = polytope(-C(j,:),-d(j)) & B;
            A1 = [temp.A -temp.b]; b1 = zeros(size(A1,1),1);

            % 2nd constraint in Equation (9) in [5]
            A2 = [eye(n), -ones(n,1); -eye(n), -ones(n,1)]; 
            b2 = zeros(2*n,1);

            % 3rd constraint in Equation (9) in [5]
            A3 = [zeros(1,n), -1]; b3 = 0;

            % 4th constraint in Equation (9) in [5]
            A4 = [eye(n), zeros(n,1)];

            % 5th constraint in Equation (9) in [5]
            A5 = [zeros(1,n), 1];

            % indices of integer variables
            intcon = [intcon size(A,2)+size(A3,2)];

            % combine constraints
            A = blkdiag(A,[A1;A2;A3]); b = [b;b1;b2;b3];
            Aeq_ = [Aeq_, [A4; A5]]; 
        end

        % combine constraints
        Aeq1 = blkdiag(Aeq1,Aeq_); Aeq2 = [Aeq2;-eye(n);zeros(1,n)];
        beq = [beq;zeros(n,1);1];
    end

    Aeq = [Aeq2,Aeq1]; A = [zeros(size(A,1),n),A];
    
    % objective function 
    f = ones(size(Aeq,2),1);

    % solve mixed-integer linear program
    options = optimoptions('intlinprog','Display','off');
    
    alpha = intlinprog(f,intcon,A,b,Aeq,beq,[],[],options);
    
    % check if a feasible solution could be found
    if isempty(alpha)
        res = true;
    else
        res = false; alpha = alpha(1:n);
    end
end

function res = aux_errorParticularSolution(eAt,U,Usum,Asum,E,dt)
% compute the error between time varying-inputs and piecewise constant
% inputs according to Proposition 2 in [3]

    res = interval(eAt*(Asum*U + E*dt*U)) + interval(eAt*Usum + E*dt*U);
end

function [Usum,Asum] = aux_precomputeAuxiliaryTerms(A,U,dt,order)
% precompute constant terms for the error between time-varying inputs and 
% piecewise constant inputs according to Proposition 2 in [3]

    % compute matrices A_i
    Apow = cell(order,1);
    Apow{1} = A*dt^2/2;

    for i = 2:order
        Apow{i} = Apow{i-1}*A*dt/(i+1);
    end

    % compute sum of all A_i in Equation (26) in [3]
    Asum = Apow{1};

    for i = 2:length(Apow)
        Asum = Asum + Apow{i};
    end

    % compute Minkowski sum of input set in Equation (26) in [3]
    Usum = Apow{1}*U;

    for i = 2:length(Apow)
        Usum = Usum + Apow{i}*U;
    end
end

function E = aux_remMatrixExp(A,dt,order)
% remainder of the matrix exponential according to Eq. (3.3) in [4]

    Aabs = abs(A);
    tmp = eye(size(A));
    Emat = expm(Aabs*dt) - tmp;
    
    for i = 1:order
       tmp = tmp * Aabs * dt / i;
       Emat = Emat - tmp;
    end
    Emat = abs(Emat);
    
    E = interval(-Emat,Emat);
end

function [F,G] = aux_curvatureEnclosure(A,E,dt,order)
% enclosure of the trajectory curvature according to Prop. 3.1 in [4]
    
    % compute auxiliary intervals
    I = cell(order+2,1);
    
    for i = 2:order+2
       I{i} = interval(i^(-i/(i-1)) - i^(-1/(i-1)),0); 
    end
    
    % compute curvature enclosure for homogeneous solution
    tmp = A*dt;
    F = E;
    
    for i = 2:order
       tmp = tmp * A * dt/i;
       F = F + I{i} * tmp;
    end

    % compute curvature enclosure for particular solution
    tmp = eye(size(A))*dt;
    G = E * dt;
    
    for i = 2:order+1
       tmp = tmp * A * dt/i;
       G = G + I{i}*tmp;
    end
end

function [T,res] = aux_constInputPropMat(A,dt)
% propagation matrix for constant inputs T = A^-1 * (e^A*dt - I)

    tmp = eye(size(A)) * dt;
    T = tmp;
    cnt = 2;

    while cnt < 1000
        tmp = tmp * dt/cnt * A;
        T = T + tmp;
        cnt = cnt + 1;
        if all(all(abs(tmp) < eps))
            break;
        end
    end

    res = cnt < 1000;
end

function list = aux_safe2unsafe(sets)
% convert a safe set defined by the union of multiple polytopes to an
% equivalent union of unsafe sets

    if ~iscell(sets)
        sets = {sets};
    end

    list = aux_reverseHalfspaceConstraints(sets{1});

    for i = 2:length(sets)

        tmp = aux_reverseHalfspaceConstraints(sets{i});

        list_ = {};

        for j = 1:length(tmp)
            for k = 1:length(list)
                if isIntersecting(list{k},tmp{j})
                    list_{end+1} = polytope(list{k}) & polytope(tmp{j});
                end
            end
        end

        list = list_;
    end
end

function res = aux_reverseHalfspaceConstraints(poly)
% get a list of reversed halfspace constraints for a given polytope

    res = {};
    poly = polytope(poly);

    for i = 1:length(poly.b)
        res{end+1} = halfspace(-poly.A(i,:),-poly.b(i));
    end
end

function pred = aux_evaluatePredicates(R,sets,index,pred)
% check if the reachable set satisfies a list of predicates defined by a 
% safe set or a union of unsafe sets 

    % loop over all reachable set/predicates combinations
    for i = 1:length(R)
        for j = 1:length(sets)
            if pred(j,i) == 2
                pred(j,i) = aux_evaluatePredicate(R{i},sets{j},index{i});
            end
        end
    end
end

function res = aux_evaluatePredicate(R,set,index)
% check if the reachable set satisfies a predicate defined by a safe set or
% a union of unsafe sets 

    if ~iscell(set)     % single safe set

        if isa(set,'polytope')
            c = set.A * R.c; 
            G = set.A * R.G;
            d = set.b;
        else
            c = set.c' * R.c; 
            G = set.c' * R.G;
            d = set.d;
        end
        G = abs(G);

        tmp = sum(G,2);

        if all(c + tmp <= d)
            res = 1;
        elseif any(c - tmp > d)
            res = 0;
        elseif any(c - sum(G(:,index),2) + ...
                                sum(G(:,setdiff(1:size(G,2),index)),2) > d)
            res = 3;
        else
            res = 2;
        end

    else                % union of unsafe sets

        list = zeros(length(set),1);

        % loop over all unsafe sets
        for i = 1:length(set)

            if isa(set{i},'polytope')
                c = set{i}.A * R.c; 
                G = set{i}.A * R.G;
                d = set{i}.b;
            else
                c = set{i}.c' * R.c; 
                G = set{i}.c' * R.G;
                d = set{i}.d;
            end

            G = abs(G);
            tmp = sum(G,2);

            if any(c - tmp > d)
                res = 1;
            elseif all(c + tmp <= d)
                res = 0; return;
            elseif all(c + sum(G(:,index),2) <= d)
                res = 3;
            else
                res = 2;

                if isa(set{i},'polytope') && ~isIntersecting(R,set{i})
                    res = 1;
                end
            end

            list(i) = res;
        end

        res = max(list(i));
    end
end

function [phi,sets] = aux_preprocessTemporalLogic(phi)
% preprocess temporal logic formula

    % convert to negation normal form
    phi = negationNormalForm(phi);

    % assign unique identifiers to all predicates
    [phi,sets] = assignIdentifiers(phi);

    % convert the regions defined by the predicates to sets
    for i = 1:length(sets)

        % convert to a union of safe sets
        eq = disjunctiveNormalForm(sets{i});
        clauses = getClauses(eq,'dnf');

        if length(clauses) == 1                 % single safe set

            sets{i} = convert2set(clauses{1});

        else                                    % union of safe sets

            list = cell(length(clauses),1);

            for j = 1:length(clauses)
                list{j} = convert2set(clauses{j});
            end

            % convert to a union of unsafe sets
            sets{i} = aux_safe2unsafe(list);
        end
    end
end

function [dt,options] = aux_initTimeStepSize(sys,options)
% choose the initial time step size

    % select a time step size that is compatible with time-varying input
    if isfield(options,'u')
        dt = options.tFinal/size(options.u,2);
    else
        dt = options.tFinal/2;
    end

    % reduce time step size until the Taylor series for the matrix
    % exponential converges
    while true
        [~,res] = aux_constInputPropMat(sys.A,dt);
        if res
            break;
        else
            dt = dt/2;
            if isfield(options,'u')
                options.u = repelem(options.u,1,2);
            end
        end
    end

    % further reduce time step size until values in the Taylor series
    % exponential remainder take reasonable values
    while true
        E = aux_remMatrixExp(sys.A,dt,10);
        if max(max(max(abs(supremum(E)))),max(max(abs(infimum(E))))) < 1e10
            break;
        else
            dt = dt/2;
            if isfield(options,'u')
                options.u = repelem(options.u,1,2);
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
