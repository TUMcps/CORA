function [A,b,Ae,be,lb,ub,intcon,ind_state,ind_rob] = robustnessEncodingMILP(phi,dt,varargin)
% robustnessEncodingMILP - encodes the robustness of a temporal logic
%   formula via liner mixed-integer constraints (see Sec. IV.D in [1])  
%
% Syntax:
%    [A,b,Ae,be,lb,ub,intcon,ind_state,ind_rob] = robustnessEncodingMILP(phi,dt)
%    [A,b,Ae,be,lb,ub,intcon,ind_state,ind_rob] = robustnessEncodingMILP(phi,dt,M)
%
% Inputs:
%    phi - temporal logic formula (class stl)
%    dt - time step for time-discretization
%    M - value for M used for the Big-M encoding for mixed-integer prog.
%
% Outputs:
%    A - matrix for the inequality constraint A*x <= b
%    b - vector for the inequality constraint A*x <= b
%    Ae - matrix for the equality constraint Ae*x = be
%    be - vector for the inequality constraint Ae*x = be
%    lb - vector specifying lower bound for the variables lb <= x <= ub
%    ub - vector specifying upper bound for the variables lb <= x <= ub
%    intcon - indizes of variables that represent integers
%    ind_state - indizes of variables that represent states of trajectory
%    ind_rob - index of the variable that represents the robustness
%
% Example: 
%    x = stl('x',2);
%    phi = finally(x(1) < 5 & x(2) < 3,interval(0,1))
%
%    [A,b,Ae,be,lb,ub,intcon,ind_state,ind_rob] = robustnessEncodingMILP(phi,0.1)
%
% References: 
%   [1] V. Raman and et al, "Model Predictive Control for Signal Temporal
%       Logic Specifications", 2016
%   [2] M. Wetzlinger and et al. "Fully Automated Verification of Linear 
%       Systems Using Inner and Outer Approximations of Reachable Sets", 
%       Transactions on Automatic Control 2023   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       24-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    M = 1e3;

    if nargin > 2
        M = varargin{1};
    end

    % initialization
    n = length(getVariables(phi));
    tFinal = maximumTime(phi); 
    N = floor(tFinal/dt)+1;
    t = (0:N)*dt;
    ind_state = reshape(1:N*n,[n,N]);

    % bring temporal logic formula to correct format and extract predicates
    [phi,~,sets] = aux_preprocessTemporalLogic(phi);

    % precompute constraints for the robustness of all predicates
    r = zeros(length(sets),N);

    len = N*n; conPred.intcon = [];
    conPred.A = zeros(1,len); conPred.b = 0; 
    conPred.Aeq = zeros(1,len); conPred.beq = 0; 
    conPred.lb = -Inf*ones(len,1); conPred.ub = Inf*ones(len,1); 

    for i = 1:length(sets)
        for j = 1:N
            [conPred,r(i,j)] = aux_constraintPredicate(conPred, ...
                                                ind_state(:,j),sets{i},M);
        end
    end

    conPred.A = conPred.A(2:end,:); 
    conPred.b = conPred.b(2:end);
    conPred.Aeq = conPred.Aeq(2:end,:); 
    conPred.beq = conPred.beq(2:end);

    % compute robustness and gradient for the STL formula
    [con,ind_rob] = aux_constraintTemporalLogic(phi,[],r,conPred,t,0,M);

    A = con.A; b = con.b; intcon = con.intcon;
    Ae = con.Aeq; be = con.beq;
    lb = con.lb; ub = con.ub;
end


% Auxiliary functions -----------------------------------------------------

function [con,r] = aux_constraintTemporalLogic(phi,con,r,conPred,time,depth,M)
% recursive function to compute the gradient of the robustness of a 
% temporal logic formula

    if ~phi.temporal

        r = r(phi.id,:);

        if isempty(con)
            con = conPred;
        end

    elseif strcmp(phi.type,'&') % ---
        % compute robustness of each hs
        [con,r1] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                conPred,time,depth+1,M);
        [con,r2] = aux_constraintTemporalLogic(phi.rhs,con,r, ...
                                                conPred,time,depth+1,M);

        r = zeros(size(r1));

        for i = 1:length(r)
            [con,r(i)] = aux_minimumConstraint([r1(i);r2(i)],M,con,1,1);
            if depth == 0
                r = r(1); break;
            end
        end

    elseif strcmp(phi.type,'|') % ---
        % compute robustness of each hs
        [con,r1] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                conPred,time,depth+1,M);
        [con,r2] = aux_constraintTemporalLogic(phi.rhs,con,r, ...
                                                conPred,time,depth+1,M);
        
        r = zeros(size(r1));

        for i = 1:length(r)
            [con,r(i)] = aux_maximumConstraint([r1(i);r2(i)],M,con,1,1);
            if depth == 0
                r = r(1); break;
            end
        end

    elseif strcmp(phi.type,'next') % ---

        [con,r] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                conPred,time,depth+1,M);

        index = find(time >= phi.from);
        index = intersect(index,1:length(r));

        r = r(index);

    elseif strcmp(phi.type,'finally') % --- 

        [con,r_] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                conPred,time,depth+1,M);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = zeros(size(r_));

        while ~isempty(index) && index(1) <= length(r)
            index = index(index <= length(r));
            [con,r(cnt)] = aux_maximumConstraint(r_(index),M,con,1,1);
            cnt = cnt + 1; index = index + 1;
            if depth == 0
                r = r(1); break;
            end
        end

    elseif strcmp(phi.type,'globally') % ---

        [con,r_] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                   conPred,time,depth+1,M);

        index = find(time >= phi.from & time <= phi.to);

        cnt = 1; 
        r = zeros(size(r_));

        while ~isempty(index) && index(1) <= length(r)
            index = index(index <= length(r));
            [con,r(cnt)] = aux_minimumConstraint(r_(index),M,con,1,1);
            cnt = cnt + 1; index = index + 1;
            if depth == 0
                r = r(1); break;
            end
        end

    elseif strcmp(phi.type,'until') % ---

        % compute robustness of each hs
        [con,r1] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                   conPred,time,depth+1,M);
        [con,r2] = aux_constraintTemporalLogic(phi.rhs,con,r, ...
                                                   conPred,time,depth+1,M);
        
        % find indices
        index = find(time >= phi.from & time <= phi.to);

        % check each index
        cnt = 1; 
        r = zeros(size(r1));
    
        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));
            r_ = zeros(1,length(index));

            % get max-min robustness
            for i = 1:length(index)
                rTmp = [r2(index(i))';r1(cnt:index(i))'];
                [con,r_(i)] = aux_minimumConstraint(rTmp,M,con,1,1);
            end

            [con,r(cnt)] = aux_maximumConstraint(r_,M,con,1,1);
    
            if depth == 0
                r = r(1); break;
            end

            cnt = cnt + 1; index = index + 1;
        end

    elseif strcmp(phi.type,'release') % ---

        % compute robustness of each hs
        [con,r1] = aux_constraintTemporalLogic(phi.lhs,con,r, ...
                                                   conPred,time,depth+1,M);
        [con,r2] = aux_constraintTemporalLogic(phi.rhs,con,r, ...
                                                   conPred,time,depth+1,M);
        
        % find indices
        index = find(time >= phi.from & time <= phi.to);
    
        % check each index
        cnt = 1; 
        r = zeros(size(r1));

        while ~isempty(index) && index(1) <= length(r)
    
            index = index(index <= length(r));

            r_ = zeros(1,length(index));
    
            for i = 1:length(index)
                % get max-min robustness
                rTmp = [r2(index(i))';r1(cnt:index(i))'];
                [con,r_(i)] = aux_minimumConstraint(rTmp,M,con,-1,1);
            end

            [con,r(cnt)] = aux_maximumConstraint(r_,M,con,1,-1);

            if depth == 0
                r = r(1); break;
            end
    
            cnt = cnt + 1; index = index + 1;
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

            safeSet = convert2set(clauses{1});
            sets{i} = aux_reverseInequalityConstraints(safeSet);

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

    % reverse first constraint
    list = aux_reverseInequalityConstraints(sets{1});

    for i = 2:length(sets)

        % reverse next constraint
        nextConstReverse = aux_reverseInequalityConstraints(sets{i});

        % go through all combinations
        list_ = {};
        for j = 1:length(nextConstReverse)
            for k = 1:length(list)
                if isa(list{k},'levelSet') || isa(nextConstReverse{j},'levelSet') || ...
                        isIntersecting_(list{k},nextConstReverse{j},'exact',1e-8)
                    % compute intersection
                    list_{end+1} = and_(list{k},nextConstReverse{j},'exact');
                end
            end
        end

        % update list
        list = list_;
    end
end

function res = aux_reverseInequalityConstraints(S)
% get a list of reversed inequality constraints for a given set

    res = {};

    if isa(S,'levelSet')
        compOp = S.compOp;

        if ~iscell(compOp)
           compOp = {compOp};
        end

        for i = 1:size(S.eq,1)
            res{end+1} = levelSet(-S.eq(i),S.vars,compOp{i});
        end

    else
        % convert to polytope
        poly = polytope(S);
        for i = 1:length(poly.b)
            res{end+1} = ~polytope(poly.A(i,:),poly.b(i));
            res{end} = normalizeConstraints(res{end},'A');
        end
    end
end

function [con,index] = aux_constraintPredicate(con,ind_x,sets,M)
% add a new variable delta with index "index" that represents the distance
% of the point x from the set S using the linear encoding in 
% Proposition 8 in [2]

    if length(sets) == 1    % single unsafe set

        S = sets{1};
        
        if size(S.A,1) == 1     % only one constraint
        
            % constraint S.A*x - S.b = delta
            con.Aeq = [con.Aeq,zeros(size(con.Aeq,1),1)];
    
            Atmp = zeros(1,size(con.Aeq,2));
            Atmp(1,ind_x) = S.A;
            Atmp(end) = -1;
            btmp = S.b;
        
            con.Aeq = [con.Aeq; Atmp];
            con.beq = [con.beq; btmp];
    
        else                    % multiple constraints
    
            % constraint S.A*x - S.b <= delta
            con.A = [con.A,zeros(size(con.A,1),1+size(S.A,1))];

            Atmp = zeros(size(S.A,1),size(con.A,2));
            Atmp(:,ind_x) = S.A;
            Atmp(:,end) = -ones(size(S.A,1),1);
            btmp = S.b;
        
            con.A = [con.A; Atmp];
            con.b = [con.b; btmp];

            % According to Eq. (3) and (5) in [1], the constraint
            %   S.A(1,:)*x - S.b(1) > delta || ... || S.A(p,:)*x - S.b(p) > delta
            % can be encoded via the following two constraint
            %   delta - S.A(i,:)*x + S.b(i) < M*(1-z(i))
            % and
            %   sum_i z(i) > 1
            % where z(i) are integer variables
            ind_z = size(con.A,2)-size(S.A,1):size(con.A,2)-1;

            Atmp = zeros(size(S.A,1),size(con.A,2));
            Atmp(:,ind_x) = -S.A;
            Atmp(:,ind_z) = M*eye(length(ind_z));
            Atmp(:,end) = ones(size(S.A,1),1);
            btmp = -S.b + M;
        
            con.A = [con.A; Atmp];
            con.b = [con.b; btmp];

            Atmp = zeros(1,size(con.A,2));
            Atmp(end,ind_z) = -ones(1,length(ind_z));
            btmp = -1;

            con.A = [con.A; Atmp];
            con.b = [con.b; btmp];

            con.ub = [con.ub;Inf*ones(1+size(S.A,1),1)];
            con.lb = [con.lb;-Inf*ones(1+size(S.A,1),1)];
            con.ub(ind_z) = 1; con.lb(ind_z) = 0;
            con.intcon = [con.intcon,ind_z];
        end

        % adapt size of the constraint matrices to the new number of variables
        con = aux_padConstraints(con);
        index = size(con.A,2);

    else            % multiple unsafe sets

        % loop over all unsafe sets
        ind_r = zeros(length(sets),1);

        for i = 1:length(sets)
            [con,ind_r(i)] = aux_constraintPredicate(con,ind_x,sets(i),M);
        end

        % take the minimum distance to any unsafe set
        [con,index] = aux_minimumConstraint(ind_r,M,con,1,1);
    end
end

function [con,index] = aux_minimumConstraint(ind,M,con,s1,s2)
% add a new variable with index "index" that represents the minimum of the 
% given variables using the MILP encoding in Eq. (3),(4) and (5) in [1]

    con.A = [con.A,zeros(size(con.A,1),1+length(ind))];
    con.Aeq = [con.Aeq,zeros(size(con.Aeq,1),1+length(ind))];
    ind_z = size(con.A,2)-length(ind):size(con.A,2)-1;
    index = size(con.A,2);

    % constraint sum_i z(i) = 1     (see Eq. (3) in [1])
    Aeq = zeros(1,size(con.Aeq,2));
    Aeq(1,ind_z) = ones(1,length(ind_z));
    
    con.Aeq = [con.Aeq;Aeq]; con.beq = [con.beq;1];
    
    % constraint x(index) <= x(ind(i))      (see Eq. (4) in [1])
    A = zeros(length(ind),size(con.A,2));
    A(:,ind) = -eye(length(ind))*s1;
    A(:,end) = ones(length(ind),1)*s2;

    con.A = [con.A;A]; con.b = [con.b;zeros(length(ind),1)];

    % constraint x(index) > x(ind(i)) - M*(1-z(i))     (see Eq. (5) in [1])
    A = zeros(length(ind),size(con.A,2));
    A(:,ind) = eye(length(ind))*s1;
    A(:,end) = -ones(length(ind),1)*s2;
    A(:,ind_z) = M*eye(length(ind));

    con.A = [con.A;A]; con.b = [con.b;M*ones(length(ind),1)];

    % constraint z(i) \in {0,1}
    con.ub = [con.ub;Inf*ones(1+length(ind),1)];
    con.lb = [con.lb;-Inf*ones(1+length(ind),1)];

    con.lb(ind_z) = 0; con.ub(ind_z) = 1;

    con.intcon = [con.intcon,ind_z];

    % adapt size of the constraint matrices to the new number of variables
    con = aux_padConstraints(con);
end

function [con,index] = aux_maximumConstraint(ind,M,con,s1,s2)
% add a new variable with index "index" that represents the maximum of the 
% given variables using the MILP encoding in Eq. (3),(4) and (5) in [1]

    con.A = [con.A,zeros(size(con.A,1),1+length(ind))];
    con.Aeq = [con.Aeq,zeros(size(con.Aeq,1),1+length(ind))];
    ind_z = size(con.A,2)-length(ind):size(con.A,2)-1;
    index = size(con.A,2);

    % constraint sum_i z(i) = 1     (see Eq. (3) in [1])
    Aeq = zeros(1,size(con.Aeq,2));
    Aeq(1,ind_z) = ones(1,length(ind_z));
    
    con.Aeq = [con.Aeq;Aeq]; con.beq = [con.beq;1];
    
    % constraint x(index) >= x(ind(i))      (see Eq. (4) in [1])
    A = zeros(length(ind),size(con.A,2));
    A(:,ind) = eye(length(ind))*s1;
    A(:,end) = -s2*ones(length(ind),1);

    con.A = [con.A;A]; con.b = [con.b;zeros(length(ind),1)];

    % constraint x(index) < x(ind(i)) + M*(1-z(i))     (see Eq. (5) in [1])
    A = zeros(length(ind),size(con.A,2));
    A(:,ind) = -eye(length(ind))*s1;
    A(:,end) = s2*ones(length(ind),1);
    A(:,ind_z) = M*eye(length(ind));

    con.A = [con.A;A]; con.b = [con.b;M*ones(length(ind),1)];

    % constraint z(i) \in {0,1}
    con.ub = [con.ub;Inf*ones(1+length(ind),1)];
    con.lb = [con.lb;-Inf*ones(1+length(ind),1)];

    con.lb(ind_z) = 0; con.ub(ind_z) = 1;

    con.intcon = [con.intcon,ind_z];

    % adapt size of the constraint matrices to the new number of variables
    con = aux_padConstraints(con);
end

function con = aux_padConstraints(con)
% adapt the constraint matrices to the new number of variables

    % determine number of variables
    len = max([size(con.A,2),size(con.Aeq,2)]);

    % adapt constraint matrices
    if size(con.A,2) < len
        con.A = [con.A,zeros(size(con.A,1),len - size(con.A,2))];
    end

    if size(con.Aeq,2) < len
        con.Aeq = [con.Aeq,zeros(size(con.Aeq,1),len - size(con.Aeq,2))];
    end

    % adapt lower and upper bound
    if length(con.lb) < len
        con.lb = [con.lb; -Inf*ones(len-length(con.lb),1)];
        con.ub = [con.ub; Inf*ones(len-length(con.ub),1)];
    end
end

% ------------------------------ END OF CODE ------------------------------
