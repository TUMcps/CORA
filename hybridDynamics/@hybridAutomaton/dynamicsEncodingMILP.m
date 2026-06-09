function [A,b,Ae,be,lb,ub,intcon,ind_state,ind_input,ind_loc,locOrig] = dynamicsEncodingMILP(HA,t,params)
% dynamicsEncodingMILP - encodes the dynamics of a hybrid automataon via 
%    linear mixed-integer constraints according to [1]  
%
% Syntax:
%    [A,b,Ae,be,lb,ub,intcon,ind_state,ind_input,ind_loc,locOrig] = dynamicsEncodingMILP(HA,t,params)
%
% Inputs:
%    HA - hybrid automaton (class hybridAutomaton)
%    t - times for the time-discretization
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
%    ind_input - indizes of variables that represent inputs of trajectory
%    ind_loc - indizes of variables representing location for trajectory
%    locOrig - assignments of fake locations to original locations
%
% Example: 
%    HA = roomHeating();
%    
%    t = 0:1:10;
%    params.R0 = interval([0;0],[3;3]);
%    params.U = interval(-5,35);
%    params.startLoc = 1;
%
%    [A,b,Ae,be,lb,ub,intcon,ind_state,ind_input] = dynamicsEncodingMILP(HA,t,params)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton
%
% References: 
%   [1] A. Bemporad and M. Morari, "Control of systems integrating logic, 
%       dynamics, and constraints", Automatica 1998

% Authors:       Niklas Kochdumper
% Written:       18-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    m = -1e5; M = 1e5; e = 1e-6; 

    % bring hybrid automaton to a suitable format
    [HA,params,locOrig] = aux_removeSelfTransitions(HA,params);
    [HA,params,locOrig_] = aux_removeDuplicateTransitions(HA,params);

    if length(locOrig_) > length(locOrig)
        for i = 1:length(locOrig)+1
            locOrig_(i) = locOrig(locOrig_(i));
        end
    end

    % initialization    
    N = length(t)-1; 

    n = HA.location(1).contDynamics.nrOfDims; 
    q = HA.location(1).contDynamics.nrOfInputs; 
    p = length(HA.location);

    % bring input arguments to the correct format
    [params,m_u,M_u] = aux_initialization(HA,params,t);

    % determine transitions with non-identity resets
    [trans,nonId] = aux_nonIdentityResets(HA);

    % initialize matrices storing indices of variables
    x = zeros(n,N+1); x(:,1) = (1:n)';
    x_ = zeros(n,N+1); x_t = zeros(n,N+1,length(nonId));
    d = zeros(p,N+1); mu = zeros(size(trans,2),N+1);
    u = zeros(q,N);
    x_m = zeros(n,N+1,p); u_m = zeros(q,N,p);

    % constraint x(t_0) \in X0
    con = setContainmentEncoding(params.R0,eye(n),zeros(n,1));

    % loop over all time steps
    for i = 1:N+1
        
        % store indices of the new variables introduced for this time step
        len = max([size(con.Ae,2),size(con.Aineq,2),length(con.lb)]);

        x_m(:,i,:) = len + reshape(1:n*p,[n,p]); len = len + n*p;
        x_t(:,i,:) = len + reshape(1:n*length(nonId),[n,length(nonId)]); len = len + n*length(nonId);
        x_(:,i) = len + (1:n)'; len = len + n;
        d(:,i) = len + (1:p)'; len = len + p;
        mu(:,i) = len + (1:size(trans,2)); len = len + size(trans,2);

        if i <= N
            u(:,i) = len + (1:q)'; len = len + q;
            u_m(:,i,:) = len + reshape(1:q*p,[q,p]); len = len + q*p;
            x(:,i+1) = len + (1:n)'; len = len + n;
        end

        % adapt size of the constraint matrices
        con.Aineq = [con.Aineq,zeros(size(con.Aineq,1),len-size(con.Aineq,2))];
        con.Ae = [con.Ae,zeros(size(con.Ae,1),len-size(con.Ae,2))]; 
        con.lb = [con.lb;-Inf*ones(len-length(con.lb),1)];
        con.ub = [con.ub;Inf*ones(len-length(con.ub),1)];
        con.intcon = [con.intcon,d(:,i)',mu(:,i)'];

        % constraint x(t_i+1) = sum_j A{j}*x_j(t_i) + B{j}*u(t_i)*d_ij + c{j}*d_ij
        if i <= N

            Ae_ = zeros(n,len); Ae_(:,x(:,i+1)) = -eye(n); 
    
            for j = 1:p
                [A,B,c] = aux_discretizedDynamics(HA.location(j), ...
                                    t(i),t(i+1),params.u{j},params.tu{j});

                Ae_(:,x_m(:,i,j)) = A;
                Ae_(:,u_m(:,i,j)) = B;
                Ae_(:,d(j,i)) = c;
            end
    
            con.Ae = [con.Ae;Ae_]; con.be = [con.be;zeros(n,1)];
        end

        % constraint x_j(t_i) = x_(t_i)*d_ij (encoded using Eq. (5b) in [1])
        for j = 1:p
            for l = 1:n
                con = aux_conProduct(con,x_m(l,i,j),x_(l,i),d(j,i),m,M);
            end
        end

        % constraint u_j(t_i) = u(t_i)*d_ij (encoded using Eq. (5b) in [1])
        if i <= N
            for j = 1:p
                for l = 1:q
                    con = aux_conProduct(con,u_m(l,i,j),u(l,i),d(j,i),m_u,M_u);
                end
            end
        end

        % constraint \sum_j d_ij = 1
        Ae_ = zeros(1,len);
        Ae_(1,d(:,i)) = ones(1,p);
        con.Ae = [con.Ae;Ae_]; con.be = [con.be;1];

        % constraint d_ij \in {0,1}
        A_ = zeros(2*p,len);
        A_(:,d(:,i)) = [eye(p);-eye(p)];
        con.Aineq = [con.Aineq;A_];
        con.bineq = [con.bineq;ones(p,1);zeros(p,1)];

        % constraint mu_ik \in {0,1}
        A_ = zeros(2*size(mu,1),len);
        A_(:,mu(:,i)) = [eye(size(mu,1));-eye(size(mu,1))];
        con.Aineq = [con.Aineq; A_];
        con.bineq = [con.bineq;ones(size(mu,1),1);zeros(size(mu,1),1)];

        % constraint P_j.A*x_j(t_i) <= d_ij*P_j.b
        for j = 1:p
            P = HA.location(j).invariant;
            A_ = zeros(size(P.A,1),len);
            A_(:,x_m(:,i,j)) = P.A; A_(:,d(j,i)) = -P.b;
            con.Aineq = [con.Aineq;A_]; 
            con.bineq = [con.bineq;zeros(size(A_,1),1)];
        end

        % constraint x_(t_i) = x(t_i) + 
        %             \sum_{k \in T} (H_k - eye(n))*xt_k(t_i) + g_k*mu_ik
        if i > 1
                
            Ae_ = zeros(n,len); 
            Ae_(:,x_(:,i)) = -eye(n);
            Ae_(:,x(:,i)) = eye(n);

            for k = 1:length(nonId)
                T = HA.location(trans(1,nonId(k))).transition(trans(2,nonId(k)));
                Ae_(:,x_t(:,i,k)) = T.reset.A - eye(n);
                Ae_(:,mu(nonId(k),i)) = T.reset.c;
            end

            con.Ae = [con.Ae;Ae_]; con.be = [con.be;zeros(n,1)];
        end

        % constraint xt_k(t_i) = x(t_i)*mu_ik (encoded using Eq. (5b) in [1])
        for k = 1:length(nonId)
            for l = 1:n
                con = aux_conProduct(con,x_t(l,i,k),x(l,i),mu(nonId(k),i),m,M);
            end
        end

        % constraint u_j(t_i) \in U_j (encoding only correct if 0 \in U)
        if i <= N
            for j = 1:p
                inputIndices = u_m(:,i,j);
                A = zeros(q,max(inputIndices)); A(:,inputIndices) = eye(q); c = zeros(q,1);
                con = setContainmentEncoding(params.U{j},A,c,con);
            end
        end

        % constraint (a_k*x(t_i) = b_k && d_i-1s(k) = 1) <=> (mu_ik = 1)
        if i > 1
            for k = 1:size(trans,2)
                T = HA.location(trans(1,k)).transition(trans(2,k));
                s = trans(1,k);
                con = aux_guardCon(con,T.guard.Ae,T.guard.be,x(:,i), ...
                                                d(s,i-1),mu(k,i),m,M,e);
            end
        end

        % constraint (d_id(k) = 1 && d_i-1s(k) = 1) <=> (mu_ik = 1)
        if i > 1
            for k = 1:size(trans,2)
                T = HA.location(trans(1,k)).transition(trans(2,k));
                s = trans(1,k);
                con = aux_conBinaryProduct(con,d(T.target,i),d(s,i-1),mu(k,i));
            end
        end

        con = padConstraints(con);
    end

    % constraint x_(t_0) = x(t_0)
    Ae_ = zeros(n,size(con.Ae,2));
    Ae_(:,x(:,1)) = eye(n);
    Ae_(:,x_(:,1)) = -eye(n);

    con.Ae = [con.Ae;Ae_]; con.be = [con.be;zeros(n,1)];

    % constraint d_1(params.startLoc) = 1
    Ae_ = zeros(1,size(con.Ae,2));
    Ae_(1,d(params.startLoc,1)) = 1;

    con.Ae = [con.Ae;Ae_]; con.be = [con.be;1];

    % assign output arguments
    A = con.Aineq; b = con.bineq; Ae = con.Ae; be = con.be;
    lb = con.lb; ub = con.ub; intcon = con.intcon;
    ind_state = x; ind_input = u; ind_loc = d;
end


% Auxiliary functions -----------------------------------------------------

function con = aux_conProduct(con,y,x,d,m,M)
% encode the constraint y = x*d with d \in {0,1} using the approach in
% Equation (5b) in [1]

    len = size(con.Aineq,2);
    A = zeros(4,len); b = zeros(4,1);

    % constraint y <= M*d
    A(1,y) = 1; A(1,d) = -M;

    % constraint y >= m*d
    A(2,y) = -1; A(2,d) = m;

    % constraint y <= x - m*(1-d)
    A(3,y) = 1; A(3,x) = -1; A(3,d) = -m; b(3) = -m;

    % constraint y >= x - M*(1-d)
    A(4,y) = -1; A(4,x) = 1; A(4,d) = M; b(4) = M;

    con.Aineq = [con.Aineq;A]; con.bineq = [con.bineq;b];
end

function con = aux_guardCon(con,a,b,x,d,mu,m,M,e)
% encode the constraint (a*x = b && d = 1) <=> (mu = 1)

    con.Ae = [con.Ae,zeros(size(con.Ae,1),4)];
 
    len = size(con.Ae,2);
    l1 = len-3; l2 = len-2; l3 = len-1; l4 = len;
    con.intcon = [con.intcon,l1,l2,l3,l4];

    % constraint l1,l2,l3,l4 \in {0,1}
    con.Aineq = blkdiag(con.Aineq,[eye(4);-eye(4)]);
    con.bineq = [con.bineq;ones(4,1);zeros(4,1)];

    % constraint (a*x <= b) <=> (l1 = 1)
    con = aux_conImplication(con,a,b,x,l1,m,M,e);

    % constraint (a*x >= b) <=> (l2 = 1)
    con = aux_conImplication(con,-a,-b,x,l2,m,M,e);

    % constraint l3 = l1*l2
    con = aux_conBinaryProduct(con,l1,l2,l3);

    % constraint l4 = l3*d
    con = aux_conBinaryProduct(con,d,l3,l4);

    % constraint l4 = mu
    Ae_ = zeros(1,size(con.Ae,2));
    Ae_(1,l4) = 1; Ae_(1,mu) = -1;

    con.Ae = [con.Ae;Ae_]; con.be = [con.be;0];
end

function con = aux_conBinaryProduct(con,d1,d2,d3)
% encode the product d3 = d1*d2 of two binary variables d1 and d2 according
% to Equation (5a) in [1]

    A = zeros(3,size(con.Aineq,2)); b = zeros(3,1);

    % constraint -d1 + d3 <= 0
    A(1,d1) = -1; A(1,d3) = 1;

    % constraint -d2 + d3 <= 0
    A(2,d2) = -1; A(2,d3) = 1;

    % constraint d1 + d2 - d3 <= 1
    A(3,d1) = 1; A(3,d2) = 1; A(3,d3) = -1; b(3) = 1;

    con.Aineq = [con.Aineq;A]; con.bineq = [con.bineq;b];
end

function con = aux_conImplication(con,a,b,x,d,m,M,e)
% encode the implication (a*x <= b) <=> (d = 1) according to Eq (4e) in [1]

    A_ = zeros(2,size(con.Aineq,2)); b_ = zeros(2,1);

    % constraint a*x - b <= M*(1-d)
    A_(1,x) = a; A_(1,d) = M; b_(1) = M+b;

    % constraint a*x - b >= e + (m-e)*d
    A_(2,x) = -a; A_(2,d) = (m-e); b_(2) = -b - e;

    con.Aineq = [con.Aineq;A_]; con.bineq = [con.bineq;b_];
end

function [params,mu,Mu] = aux_initialization(HA,params,t)
% initialize some useful auxiliary variables

    % create zero input set in case input set is not provided
    if ~isfield(params,'U')
        params.U = cell(length(HA.location),1);
        for i = 1:length(HA.location)
            m = HA.location(i).contDynamics.nrOfInputs;
            params.U{i} = zonotope(zeros(m,1));
        end
    else
        if ~iscell(params.U)
            params.U = repmat({params.U},[length(HA.location),1]);
        end
    end

    % consider time-varying inputs
    if ~isfield(params,'u')
        params.u = cell(length(HA.location),1);
        for i = 1:length(HA.location)
            m = HA.location(i).contDynamics.nrOfInputs;
            params.u{i} = zeros(m,1);
        end
    else
        if ~iscell(params.u)
            params.u = repmat({params.u},[length(HA.location),1]);
        end
    end

    % generate time points for time-varying inputs
    params.tu = cell(length(HA.location),1);

    for i = 1:length(HA.location)
        params.tu{i} = linspace(t(1),t(end),size(params.u{i},2)+1);
    end

    % center input sets at zero, such that 0 \in U
    for i = 1:length(HA.location)
        c = center(params.U{i});
        params.u{i} = params.u{i} + c;
        params.U{i} = params.U{i} + (-c);
    end

    % get upper and lower bound on the input such that mu <= u <= Mu
    Mu = 0;

    for i = 1:length(params.U)
        I = interval(params.U{i});
        inputBounds = [abs(infimum(I));abs(supremum(I))];
        Mu = max(Mu,max(inputBounds));
    end

    mu = -Mu;
end

function [A,B,c] = aux_discretizedDynamics(loc,tStart,tEnd,u,tu)
% compute the time-discretized dynamcis

    % find corresponding inputs
    t = sort(unique([tStart,tu(tu >= tStart & tu <= tEnd),tEnd]));
    u = interp1(tu',[u,u(:,end)]',t','previous')';

    % loop over all time steps
    n = loc.contDynamics.nrOfDims; m = loc.contDynamics.nrOfInputs;
    A = eye(n); B = zeros(n,m); c = zeros(n,1);

    for i = 1:length(t)-1
        dt = t(i+1) - t(i);
        sys = linearSysDT(loc.contDynamics,dt);
        A = sys.A*A;
        B = sys.A*B + sys.B;
        c = sys.A*c + sys.B*u(:,i) + sys.c;
    end

end

function [trans,nonId] = aux_nonIdentityResets(HA)
% determine all transitions that have non-identity reset functions

    trans = []; nonId = [];

    % loop over all transitions
    for i = 1:length(HA.location)
        for j = 1:length(HA.location(i).transition)

            T = HA.location(i).transition(j);
            n = length(T.reset.c);
            trans = [trans,[i;j]];

            if sum(sum(abs(T.reset.A - eye(n)))) > eps || sum(abs(T.reset.c)) > eps
                nonId = [nonId,size(trans,2)];
            end
        end
    end
end

function [HA,params,locOrig] = aux_removeSelfTransitions(HA,params)
% remove all transitions to the same mode from the hybrid automaton

    % find first self transition in the hybrid automaton
    found = false; locOrig = 1:length(HA.location);

    for i = 1:length(HA.location)
        trans = HA.location(i).transition;
        for j = 1:length(trans)
            if trans(j).target == i
                found = true; ind = [i,j]; locOrig = [locOrig,i]; break;
            end
        end

        if found
            break;
        end
    end

    % no self transition -> finished
    if ~found
        return;

    % remove self transition by adding additional state
    else
        locs = HA.location; L = locs(ind(1));
        trans = L.transition; T = trans(ind(2));

        % introduce additional "fake" location
        locs = [locs;location(L.invariant,trans,L.contDynamics)];

        % modify original location
        trans(ind(2)) = transition(T.guard,T.reset,length(locs));
        locs(ind(1)) = location(L.invariant,trans,L.contDynamics);
        
        HA = hybridAutomaton(locs);

        % modify parameters
        if isfield(params,'u') && iscell(params.u)
            params.u{end+1} = params.u{ind(1)};
        end

        if isfield(params,'U') && iscell(params.U)
            params.U{end+1} = params.U{ind(1)};
        end
    end

    % call function again until no self transitions are left
    [HA,params,locOrig_] = aux_removeSelfTransitions(HA,params);

    if length(locOrig_) > length(locOrig)
        for i = 1:length(locOrig)+1
            locOrig_(i) = locOrig(locOrig_(i));
        end
    end
end

function [HA,params,locOrig] = aux_removeDuplicateTransitions(HA,params)
% remove all duplicate transitions that start and end at the same source
% and destination modes

    % find first duplicate transition in the hybrid automaton
    found = false; locOrig = 1:length(HA.location);

    for i = 1:length(HA.location)
        trans = HA.location(i).transition;
        for j = 1:length(trans)
            for k = j+1:length(trans)
                if trans(j).target == trans(k).target
                    found = true; ind = [i,k]; locOrig = [locOrig,i];
                    break;
                end
            end
        end

        if found
            break;
        end
    end

    % no duplicate transitions -> finished
    if ~found
        return;

    % remove duplicate transition by adding additional state
    else
        locs = HA.location; L = locs(ind(1));
        trans = L.transition; T = trans(ind(2)); 
        L2 = locs(T.target); trans2 = L2.transitions;

        % introduce additional "fake" location
        locs = [locs;location(L2.invariant,trans2,L2.contDynamics)];

        % modify original location
        trans(ind(2)) = transition(T.guard,T.reset,length(locs));
        locs(ind(1)) = location(L.invariant,trans,L.contDynamics);
        
        HA = hybridAutomaton(locs);

        % modify parameters
        if isfield(params,'u') && iscell(params.u)
            params.u{end+1} = params.u{ind(1)};
        end

        if isfield(params,'U') && iscell(params.U)
            params.U{end+1} = params.U{ind(1)};
        end
    end

    % call function again until no self transitions are left
    [HA,params,locOrig_] = aux_removeDuplicateTransitions(HA,params);

    if length(locOrig_) > length(locOrig)
        for i = 1:length(locOrig)+1
            locOrig_(i) = locOrig(locOrig_(i));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
