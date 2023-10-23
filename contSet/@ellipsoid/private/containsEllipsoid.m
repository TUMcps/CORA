function res = containsEllipsoid(E1,E2)
% containsEllipsoid - checks whether an ellipsoid contains another ellipsoid
%
% Syntax:
%    E = containsEllipsoid(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object
%
% Outputs:
%    res - true/false
%
% References:
%    [1] Boyd et al. Convex Optimization (B.2, ex. B.1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       16-March-2021
% Last update:   27-July-2021 (allowed both E1,E2 to be degenerate)
%                23-May-2022 (solver-specific implementation of SDP to avoid yalmip)
%                06-July-2022 (VG, rescaling of M/Ms to avoid numerical issues)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = dim(E1);

% normalize to prevent numerical issues:
max_val = max([svd(E1.Q);svd(E2.Q)]);
E1 = 1/sqrt(max_val)*E1;
E2 = 1/sqrt(max_val)*E2;
q1 = E1.q;
E1 = E1 + (-q1);
E2 = E2 + (-q1);

% if both are degenerate
if ~isFullDim(E1) && ~isFullDim(E2)
    nt1 = rank(E1);
    nt2 = rank(E2);
    % cannot be contained if one has more actual dimensions
    if nt1<nt2
        res = false;
        return;
    end
    Q_sum = E1.Q+E2.Q;
    n_nd = rank(ellipsoid(Q_sum));
    % can only be contained if degeneracies are aligned
    if n_nd<n
        [T,~,~] = svd(Q_sum);
        % axis-align degeneracies
        E1t = T'*E1;
        E2t = T'*E2;
        if all(withinTol(E1t.q(n_nd+1:end),E2t.q(n_nd+1:end),E1t.TOL))
            res = containsEllipsoid(project(E1t,1:n_nd),project(E2t,1:n_nd));
        else
            res = false;
        end
        return;
    else
        res = false;
        return;
    end
    
elseif ~isFullDim(E1) && isFullDim(E2)% E1 is the degenerate one, E2 cannot be contained
    res = false;
    return;
elseif isFullDim(E1) && ~isFullDim(E2) % E2 is degenerate
    nt = rank(E2);
    if nt==0
        % E2 contains only single point
        res = contains_(E1,E2.q,'exact',0);
        return;
    end
    [T,~,~] = svd(E2.Q);
    E2 = T'*E2;
    E1 = T'*E1;
    % project
    x2_rem = E2.q(nt+1:end);
    E2 = project(E2,1:nt);
    % resolve x2_rem in E1t by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    n = dim(E1);
    I = eye(n);
    for i=1:(n-nt)
        Hi = conHyperplane(I(nt+i,:),x2_rem(i));
        E1 = and_(E1,Hi,'outer');
        % if intersection is empty, E2 cannot be contained in E1
        if representsa_(E1,'emptySet',eps)
            res = false;
            return;
        end
    end
    % now, E1 also has "0"s at nt+1:end
    E1 = project(E1,1:nt);
end

% E1, E2 non-degenerate

% both full-dimensional
TOL = min(E1.TOL,E2.TOL);

if all(withinTol(E1.q,E2.q,TOL)) 
    res = isBigger(E1,E2);
    return;
end

Qs = E1.Q;
Q = E2.Q;
qs = E1.q;
q = E2.q;

Q_inv = inv(Q);
Qs_inv = inv(Qs);
M = [Q_inv, -Q_inv*q;-(Q_inv*q)', q'*Q_inv*q-1];
Ms = [Qs_inv, -Qs_inv*qs;-(Qs_inv*qs)', qs'*Qs_inv*qs-1];

% rescale
fac = max(abs(M),[],'all');
M = 1/fac*M;
Ms = 1/fac*Ms;

n = size(Qs,1); % Q and Qs have same dimensions (see above)

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end
persistent isSDPT3
if isempty(isSDPT3)
    isSDPT3 = isSolverInstalled('sdpt3');
end


if isMosek
    % model the problem ourselves using MOSEK
    % see @ellipsoid/private/andEllipsoidIA.m for general MOSEK problem
    % description

    % For this specific problem, we identify these matrices using the dual
    % problem formulation (easier)

    % in our case, we want to solve the problem
    %%% t*M >= M_s,
    % where M = [Q_inv, -Q_inv*q;-(Q_inv*q)', q'*Q_inv*q-1] and
    %  Ms = [Qs_inv, -Qs_inv*qs;-(Qs_inv*qs)', qs'*Qs_inv*qs-1];
    % Thus, we only need to check feasibility, and the problem is already
    % given in dual form, i.e.
    %%% find t 
    %%%     s.t.    -M_s -(-M*t) >= 0,
    %%%                        t >= epsilon.
    
    % Since dim(M) = n+1, so is the dimension of the single PSD matrix
    prob.bardim = n+1;
    
    [rl_ms,cl_ms,val_ms] = find(tril(-Ms));
    [rl_m,cl_m,val_m] = find(tril(-M));
    
    N = (n+1)*(n+1+1)/2;
    % Cb = -M_s
    prob.barc.subj = ones(1,length(rl_ms));
    prob.barc.subk = rl_ms';
    prob.barc.subl = cl_ms';
    prob.barc.val = val_ms';

    % Ab = M
    prob.bara.subi = ones(1,length(rl_m));
    prob.bara.subj = ones(1,length(rl_m));
    prob.bara.subk = rl_m';
    prob.bara.subl = cl_m';
    prob.bara.val = val_m';
    
    % we have exactly one constraint, and zero scalar/cone variables
    % => size(prob.a) = [1,0]
    % scalar constraints part
    % but we do have one lower-bound constraint
    prob.a = sparse(zeros(1,0),zeros(1,0),zeros(1,0),1,0);

    % objective in scalar variables
    prob.c = [];

    % inequality constraints (equality constraints just have same lb and ub)
    prob.blc = 0;
    prob.buc = 0;

    % optimize
    [~,res] = mosekopt('minimize echo(0)',prob);

    % check if everything worked out
    if res.rcode==0 
        % all good
        % check if feasible
        if ~strcmp(res.sol.itr.solsta,'OPTIMAL')
            % either infeasible (which is fine), or some funny business
            % if solsta and prosta contain "infeasible" (little ugly),
            % assume infeasibility
            if contains(res.sol.itr.prosta,'INFEASIBLE') && ...
               contains(res.sol.itr.solsta,'INFEASIBLE')
                res = false;
            else
                % this might be a little too harsh
                throw(CORAerror('CORA:solverIssue','mosek'));
            end
        else
            res = true;
        end
    else
        throw(CORAerror('CORA:solverIssue','mosek'));
    end

elseif isSDPT3
    % IMPORTANT: Vectorization is upper-triangular (in contrast to MOSEK,
    % which is lower-triangular)
    % ALSO: For some reason, off-diagonal elements when stacking upper- 
    % triangular matrices should be multiplied by sqrt(2)...
    blk = {'s',1+n};
    
    % coefficient matrix C (mosek: Cb)
    C{1} = sparse(-Ms);
    
    % y = t
    n_y = 1;

    nM = (n+1)*(n+1+1)/2;

    % constraint matrix: -M
    % multiply off-diagonal entries by sqrt(2)
    M_ = sqrt(2)*(M-eye(n+1).*M) + eye(n+1).*M;
    % get index for upper triangular matrix
    a_subk = 1:nM;
    a_subl = ones(1,nM);
    a_val = -M_(triu(ones(n+1))==1)';
    
    At{1} = sparse(a_subk,a_subl,a_val,nM,n_y);

    % since sdpt3 supports logdet explicitly, b = 0
    b = 0;
    
    OPTIONS.printlevel = 0;

    % call solver
    [~,~,~,~,info] = sqlp(blk,At,C,b,OPTIONS);
    
    if info.termcode~=0
        % either primal and/or dual infeasible (should only ever occur
        % both), or numerical issues
        if any(info.termcode==[1,2])
            res = false;
        else
            % do not know error code, assuming problem
            throw(CORAerror('CORA:solverIssue','sdpt3'));
        end
    else
        res = true;
    end
    
elseif isYalmipInstalled()
    %For more details, see [1]
    t = sdpvar;
    
    options = sdpsettings('verbose',0);
    %when solving problems with higher-dimensional ellipsoids, sdpt3 [2] is
    %recommended to avoid numerical issues
    if exist('sdpt3','file')
        options.solver = 'sdpt3';
    end
    diagnostics = optimize([t*M>=Ms,t>=0],[],options);
    %either feasible or not feasible
    if ~any(diagnostics.problem==[1,0])
         throw(CORAerror('CORA:solverIssue'));
    end
    if value(t)==0
        throw(CORAerror('CORA:specialError',...
            'problem with solution of feasibility problem (t)'));
    end
    res = ~diagnostics.problem;
    warning("YALMIP was used to model the problem - " + ...
        "consider installing a supported solver to speed up computation...");

else
    throw(CORAerror('CORA:noSuitableSolver','SDP'));
end

% ------------------------------ END OF CODE ------------------------------
