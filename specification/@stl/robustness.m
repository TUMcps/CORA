function val = robustness(obj,varargin)
% robustness - compute robustness of a temporal logic formula according to
%              Definition 2.3 in [1]
% 
% Syntax:
%    val = robustness(obj,R)
%    val = robustness(obj,sim)
%    val = robustness(obj,x,t)
%
% Inputs:
%    obj - logic formula (class stl)
%    R - reachable set (class reachSet)
%    sim - simulated traces (class simResult)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%
% Outputs:
%    val - robustness of the STL formula (val > 0 if satisfied, val < 0 if
%          violated)
%
% Example: 
%    x = stl('x',2);
%    eq = until(x(2) < -0.5,x(1) > 0.5,interval(0,1));
%    
%    phi = -pi/2:0.01:0;
%    x = [cos(phi'),sin(phi')];
%    t = linspace(0,1,length(phi))';
%
%    val = robustness(eq,x,t)
%
% References:
%    [1] Sankaranarayanan, Sriram., et al. "Falsification of Temporal 
%        Properties of Hybrid Systems Using the Cross-Entropy Method", 
%        HSCC 2012.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl/modelCheckTrace

% Authors:       Niklas Kochdumper
% Written:       23-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    if nargin > 2
        x = varargin{1};
        t = varargin{2};
        type = 'trace';
    elseif isa(varargin{1},'reachSet')
        R = varargin{1};
        type = 'reachSet';
    elseif isa(varargin{1},'simResult')
        sim = varargin{1};
        type = 'simResult';
    else
        throw(CORAerror('CORA:wrongValue', "second", ...
                 'must by of class "reachSet" or "simResult".'));
    end

    % bring temporal logic formula to correct format and extract predicates
    [phi,pred,sets] = aux_preprocessTemporalLogic(obj);

    % precompute the robustness for all predicates
    if strcmp(type,'trace')
        [r_pred,t] = aux_robustnessTrace(x,t,sets);
    elseif strcmp(type,'reachSet')
        [r_pred,t,phi] = aux_robustnessReachSet(R,phi,sets,pred);
    else
        val = aux_robustnessSimResult(sim,phi,sets);
        return;
    end
    
    % compute robustness of the overall temporal logic formula
    r = aux_robustnessTemporalLogic(phi,r_pred,t);
    val = r(1);
end


% Auxiliary functions -----------------------------------------------------

function r = aux_robustnessTemporalLogic(phi,r,time)
% recursive function to compute the robustness of a temporal logic formula
% according to Definition 2.3 in [1]

    if ~phi.temporal

        r = r(phi.id,:);

    elseif strcmp(phi.type,'&')

        r1 = aux_robustnessTemporalLogic(phi.lhs,r,time);
        r2 = aux_robustnessTemporalLogic(phi.rhs,r,time);

        r = min([r1;r2],[],1);

    elseif strcmp(phi.type,'|')

        r1 = aux_robustnessTemporalLogic(phi.lhs,r,time);
        r2 = aux_robustnessTemporalLogic(phi.rhs,r,time);

        r = max([r1;r2],[],1);

    elseif strcmp(phi.type,'next')

        r = aux_robustnessTemporalLogic(phi.lhs,r,time);

        index = find(time >= phi.from);
        index = intersect(index,1:length(r));

        r = [r(index), -inf*ones(1,length(r)-length(index))];

    elseif strcmp(phi.type,'finally')

        r_ = aux_robustnessTemporalLogic(phi.lhs,r,time);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r_));

        while ~isempty(index) && index(1) <= length(r)
            index = index(index <= length(r));
            r(cnt) = max(r_(index));
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'globally')

        r_ = aux_robustnessTemporalLogic(phi.lhs,r,time);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r_));

        while ~isempty(index) && index(1) <= length(r) 
            index = index(index <= length(r));
            r(cnt) = min(r_(index));
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'until')

        r1 = aux_robustnessTemporalLogic(phi.lhs,r,time);
        r2 = aux_robustnessTemporalLogic(phi.rhs,r,time);
        
        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r1));
    
        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));
    
            for i = 1:length(index)
                r(cnt) = max(r(cnt),min(r2(index(i)),min(r1(cnt:index(i)))));
            end
    
            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'release')

        r1 = aux_robustnessTemporalLogic(phi.lhs,r,time);
        r2 = aux_robustnessTemporalLogic(phi.rhs,r,time);
        
        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = -inf * ones(size(r1));
        r1 = -r1; r2 = -r2;
    
        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));
    
            for i = 1:length(index)
                r(cnt) = max(r(cnt),min(r2(index(i)),min(r1(cnt:index(i)))));
            end
    
            cnt = cnt + 1; index = index + 1;
        end

        r(~isinf(r)) = -r(~isinf(r));
    end
end

function [r,t] = aux_robustnessTrace(x,t,sets)
% precompute the robustness for all predicates on a single trace

    % bring to common time step size
    dt = min(diff(t));

    if abs(dt - max(diff(t))) > eps
        t_ = t;
        t = 0:dt:t(end);
        [~,ind] = unique(t_);
        x = [interp1(t_(ind),x(ind,1),t,'linear','extrap'); ...
             interp1(t_(ind),x(ind,2),t,'linear','extrap')]';
    end

    % compute robustness
    r = zeros(length(sets),length(t));

    for i = 1:length(t)
        for j = 1:length(sets)
            r(j,i) = aux_robustnessPoint(x(i,:)',sets{j});
        end
    end
end

function [r,t,eq] = aux_robustnessReachSet(R,eq,sets,pred)
% precompute the robustness for all predicates on the reachable sets

    % determine time step size
    [dt,uniform,hybrid] = timeStepSize(R);

    % construct time vector
    tStart = R(1).timePoint.time{1};

    if isa(tStart,'interval')
        tStart = infimum(tStart);
    end

    tFinal = min(maximumTime(eq),R.timePoint.time{end});

    if ~uniform
        t = linspace(tStart,tFinal,ceil((tFinal-tStart)/min(dt)));
        dt = (tFinal-tStart)/(length(t)-1);
    else
        t = tStart:dt:max(tFinal,dt);

        if t(end) < tFinal - eps
            t = [t, t(end)+dt];
        end
    end

    % make temporal logic formula consistent with time step size
    eq = consistentTimeStep(eq,dt);

    % compute robustness for all predicates on the reachable set
    r = zeros(length(sets),length(t)-1);

    for j = 1:length(sets)

        % check if predicate is identical to another predicate
        exists = false;

        for i = 1:j-1
            if isequal(pred{i},pred{j})
                r(j,:) = r(i,:);
                exists = true; break;
            end
        end

        if exists
            continue;
        end

        % compute robustness
        for i = 1:length(t)-1

            if ~uniform || hybrid || size(R,1) > 1
                R_ = aux_findReachSets(R,interval(t(i),t(i+1)));
            else
                R_ = R.timeInterval.set(i);
            end

            for k = 1:length(R_)
                r(j,i) = aux_robustnessSet(R_{k},sets{j});
                if ~r(j,i)
                    break;
                end
            end
        end
    end
end

function val = aux_robustnessSimResult(sim,phi,sets)
% compute robustness for a set of traces (= minimum over all traces)

    val = inf;

    for i = 1:length(sim)
        for j = 1:length(sim(i,1))
            [r_pred,t] = aux_robustnessTrace(sim(i,1).x{j},sim(i,1).t{j},sets);
            r = aux_robustnessTemporalLogic(phi,r_pred,t);
            val = min(val,r(1));
        end
    end
end

function val = aux_robustnessPoint(p,set)
% compute robustness of a single point (=min. distance from the unsafe set)

    if ~iscell(set)     % single unsafe set

        if isa(set,'polytope')
            if contains(set,p)
                len = sqrt(sum(set.A.^2,2));
                val = max((set.A * p - set.b)./len);
            else
                n = length(p);
                options = optimoptions('quadprog','display','off');
                p_ = quadprog(eye(n),-p,set.A,set.b, ...
                                                [],[],[],[],[],options);
                val = norm(p-p_);
            end
        elseif isa(set,'halfspace')
            val = (set.c'*p - set.d)/norm(set.c);
        else
            val = max(set.funHan(p));
        end

    else                % union of unsafe sets

        val = inf;

        % loop over all unsafe sets
        for i = 1:length(set)
            val = min(val,aux_robustnessPoint(p,set{i}));
        end
    end
end

function val = aux_robustnessSet(S,set)
% compute robustness of a single point (=min. distance from the unsafe set)

    if ~iscell(set)     % single unsafe set

        if isa(set,'polytope')
            S = zonotope(S);
            if ~isIntersecting_(set,S,'approx') || ...
                                    ~isIntersecting_(set,S,'exact')

                % solve quadratic program with variables [d;x;\alpha]
                %
                %   min ||d||^2 s.t. d = c + G*\alpha - x,
                %                    A*x <= b,
                %                    -1 <= \alpha <= 1
        
                A = set.A; b = set.b; n = size(A,2);
                c = center(S); G = generators(S); m = size(G,2);
        
                % constraint d = c + G*\alpha - x
                Ceq = [eye(n),eye(n),-G]; deq = c;
        
                % constraint A*x <= b
                C1 = [zeros(size(A,1),n),A]; d1 = b;
        
                % constraint -1 <= \alpha <= 1
                C2 = [eye(m);-eye(m)]; d2 = ones(2*m,1);
        
                % combined inequality constraints
                C = blkdiag(C1,C2); d = [d1;d2];
        
                % objective function
                H = blkdiag(2*eye(n),zeros(n+m));
                f = zeros(2*n+m,1);
        
                % solve quadratic program
                options = optimoptions('quadprog','display','off');
        
                [~,val] = quadprog(H,f,C,d,Ceq,deq,[],[],[],options);
        
                val = sqrt(val);
                
            else
                % solve linear program with variables [x;r;\alpha]
                %
                %   max r s.t. A(i,:)*x + ||A(i,:)||*r <= b,
                %              x = c + G * \alpha,
                %              r >= 0,
                %              -1 <= \alpha <= 1
        
                A = set.A; b = set.b; n = size(A,2);
                c = center(S); G = generators(S); m = size(G,2);
        
                % constraint A(i,:)*x+||A(i,:)||*r <= b
                C1 = [A,sqrt(sum(A.^2,2))]; d1 = b;
        
                % constraint r >= 0
                C2 = [zeros(1,size(A,2)),-1]; d2 = 0;
        
                % constraint -1 <= \alpha <= 1
                C3 = [eye(m);-eye(m)]; d3 = ones(2*m,1);
        
                % constraint x = c + G * \alpha,
                Ceq = [eye(n),zeros(n,1),-G]; deq = c;
        
                % combined inequality constraints
                C = blkdiag([C1;C2],C3); d = [d1;d2;d3];
        
                % objective function
                f = zeros(size(Ceq,2),1); f(n+1) = -1;
        
                % solve linear program
                options = optimoptions('linprog','display','off');
        
                [~,val] = linprog(f,C,d,Ceq,deq,[],[],options);
            end

        elseif isa(set,'halfspace')
            val = (supportFunc(S,set.c,'lower') - set.d)/norm(set.c);
        else
            val = max(infimum(set.funHan(interval(S))));
        end

    else                % union of unsafe sets

        val = inf;

        % loop over all unsafe sets
        for i = 1:length(set)
            val = min(val,aux_robustnessSet(S,set{i}));
        end
    end
end

function [phi,pred,sets] = aux_preprocessTemporalLogic(phi)
% preprocess temporal logic formula

    % convert to negation normal form
    phi = negationNormalForm(phi);

    % assign unique identifiers to all predicates
    [phi,pred] = assignIdentifiers(phi);

    % convert the regions defined by the predicates to sets
    sets = cell(size(pred));

    for i = 1:length(pred)

        % convert to a union of safe sets
        eq = disjunctiveNormalForm(pred{i});
        clauses = getClauses(eq,'dnf');

        if length(clauses) == 1                 % single safe set

            tmp = convert2set(clauses{1});
            sets{i} = aux_reverseInequalityConstraints(tmp);

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

function list = aux_safe2unsafe(sets)
% convert a safe set defined by the union of multiple sets to an
% equivalent union of unsafe sets

    list = aux_reverseInequalityConstraints(sets{1});

    for i = 2:length(sets)

        tmp = aux_reverseInequalityConstraints(sets{i});

        list_ = {};

        for j = 1:length(tmp)
            for k = 1:length(list)
                if isa(list{k},'levelSet') || isa(tmp{j},'levelSet') || ...
                        isIntersecting_(list{k},tmp{j},'exact')
                    if isa(list{k},'halfspace') && isa(tmp{j},'halfspace')
                        list_{end+1} = and_(polytope(list{k}),...
                                                    polytope(tmp{j}),'exact');
                    else
                        list_{end+1} = and_(list{k},tmp{j},'exact');
                    end
                end
            end
        end

        list = list_;
    end
end

function res = aux_reverseInequalityConstraints(set)
% get a list of reversed inequality constraints for a given set

    res = {};

    if isa(set,'levelSet')

        compOp = set.compOp;

        if ~iscell(compOp)
           compOp = {compOp};
        end

        for i = 1:size(set.eq,1)
            res{end+1} = levelSet(-set.eq(i),set.vars,compOp{i});
        end

    elseif isa(set,'halfspace')

        res = {halfspace(-set.c,-set.d)};

    else

        poly = polytope(set);
    
        for i = 1:length(poly.b)
            res{end+1} = polytope(-poly.A(i,:),-poly.b(i));
        end
    end
end

function sets = aux_findReachSets(R,time)
% get all sets that belong to the given time

    sets = {};

    % get all reachable sets that belong to the given time
    R = find(R,'time',time);

    % loop over all reachable set objects
    for i = 1:size(R,1)

        if isnumeric(time) && ~isempty(R(i).timePoint)
            sets = [sets; R(i).timePoint.set];
        else
            if ~isempty(R(i).timePoint)
                sets = [sets; R(i).timePoint.set];
            end
            if ~isempty(R(i).timeInterval)
                sets = [sets; R(i).timeInterval.set];
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
