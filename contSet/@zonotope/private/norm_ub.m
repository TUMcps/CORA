function ub = norm_ub(Z,type)
% norm_ub - computes bound on the maximum norm value
%
% Syntax:
%    ub = norm_ub(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    ub - upper bound on the maximum norm value
%
% Example: 
%    Z = zonotope([0;0],[1 3 -2; 2 -1 0]);
%    norm(Z,2,'ub_convex')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: norm

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   23-May-2022 (VG, model optimization problem directly)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~exist('type','var')
    type = 2;
end
if type~=2
    throw(CORAerror('CORA:notSupported','Only Euclidean norm supported'));
end
if ~all(Z.c==0)
    throw(CORAerror('CORA:notSupported','Not implemented for non-zero center'));
end

G = Z.G;
[~,m] = size(G);
GG = G'*G;

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end
persistent isSDPT3
if isempty(isSDPT3)
    isSDPT3 = isSolverInstalled('sdpt3');
end

%compute upper bound on norm via dual problem of max_{|u|<=1} u'*G'*G*u
if isMosek
    % The dual problem to 
    %%% max_{|u|<=1} u'*G'*G*u
    % is given by
    %%% min_{d>=0}  ones(1,m)*d,
    %%%     s.t.   diag(d) - G'*G >= 0. 

    % The (reduced; not everything need) standard primal SDP form using
    % MOSEK is (A = [--a^{(i),T}--]; At = [...,a^{(i)},...])
    %%% min_{X,x}     C#X,
    %%%     s.t.    a^{(i),T}*x + Ai#X <= ub_c_i,
    %%%             -a^{(i),T}*x - Ai#X <= -lb_c_i,
    %%%             x <= ub_x,
    %%%             -x <= -lb_x,
    %%%             X >= 0.
    % The corresponding dual problem is given by
    %%% max_{mxl,mxu,mcl,mcu}   -ub_c'*mcu + lb_c'*mcl - ub_x'*mxu +
    %%%                         lb_x'*mxl,
    %%%  (1)     s.t.            At*(mcu-mcl) + mxu-mxl = 0,
    %%%                         C + sum_{i=1}^{m} Ai*(mcu-mcl)_i >= 0,
    %%%                         mxl,mxu,mcl,mcu >= 0.
    %%%
    % Transferred to the dual domain, our problem to solve is given by
    %%% min_{d}     ones(1,m)*d,
    %%%     s.t.    -G'*G + sum_{i=1}{m}e_i*e_i'*d_i >= 0,
    %%%             d >= 0.
    % Now, set ub_c = lb_c, ub_x = inf(m,1), and let lb_x = zeros(m,1). 
    % Further, set At = eye(m). ub_x forces mxu = 0,
    % and the choice of lb_x makes the choice of mxl irrelevant for the
    % objective value. Lastly, let ub_c=lb_c = ones(m,1). Thus, (1) becomes
    %%%-min_{d}     ones(m,1)'*d,
    %%%     s.t.    -G'*G + sum_{i=1}{m}e_i*e_i'*d_i >= 0,
    %%%             d = mxl >= 0,
    % where d = mcu-mlc.
    
    prob.bardim = m;

    [rl_GG,cl_GG,val_GG] = find(tril(-GG));
    [rl_A,cl_A,val_A] = find(eye(m));
    
    % Cb = -G'*G
    prob.barc.subj = ones(1,length(val_GG));
    prob.barc.subk = rl_GG';
    prob.barc.subl = cl_GG';
    prob.barc.val = val_GG';

    % Ab_i = -e_i*e_i'
    prob.bara.subi = 1:m;
    prob.bara.subj = ones(1,m);
    prob.bara.subk = rl_A';
    prob.bara.subl = cl_A';
    prob.bara.val = val_A';

    % scalar constraints part: no scalar variables
    prob.a = sparse(eye(m));

    % objective in scalar variables
    prob.c = [];

    % inequality constraints (equality constraints just have same lb and ub)
    prob.blc = ones(1,m);
    prob.buc = prob.blc;

    % d>=0
    prob.blx = zeros(1,m);
    prob.bux = inf(1,m);

    % optimize
    [~,res] = mosekopt('minimize echo(0)',prob);

    % check if everything worked out
    if res.rcode~=0 || ~strcmp(res.sol.itr.solsta,'OPTIMAL')
        % this might be a little too harsh
        throw(CORAerror('CORA:solverIssue','mosek'));
    end
    
    % extract solution (d = mcu-mcl)
    d_sol = res.sol.itr.suc-res.sol.itr.slc;
    ub = sqrt(d_sol'*ones(m,1));


elseif isSDPT3
    % solve using sdpt3
    blk = cell(2,2);
    % first block: linear cone (more or less x>=0)
    blk(1,:) = {'l',m};
    % second block sdp
    blk(2,:) = {'s',m};

    % coefficient matrix C (mosek: Cb)
    C{2} = sparse(-GG);

    % block for d >= 0
    C{1} = sparse([],[],[],m,1);

    % At (sdp part): -G'*G - sum -e_i*e_i'*d_i >= 0
    gsum = @(x) x.*(x+1)./2;
    asubk = gsum(1:m);
    asubl = 1:m;
    aval = -ones(1,m);
    At{2} = sparse(asubk,asubl,aval,gsum(m),m);

    % At (linear cone part): At{1}*y + z_l = 0 (C{1} = 0)
    % Thus, with At{1} = -eye(m) => d = z_l >= 0.
    At{1} = sparse(-eye(m));
    
    % since dual has max => -ones, so we have ~ min d'*ones(m,1)
    b = -ones(m,1);
    
    OPTIONS.printlevel = 0;

    % call solver
    [~,~,y_sol,~,info] = sqlp(blk,At,C,b,OPTIONS);

    if info.termcode~=0
        throw(CORAerror('CORA:solverIssue','sdpt3'));
    end

    % extract solution
    d_sol = y_sol;
    ub = sqrt(d_sol'*ones(m,1));

elseif isYalmipInstalled()
    d = sdpvar(m,1);
    options = sdpsettings('verbose',0);
    %compute optimization problem
    optimize([GG<=diag(d),d>=0],d'*ones(m,1),options);
    ub = sqrt(value(d'*ones(m,1)));
    warning("YALMIP was used to model the problem - consider installing a supported solver to speed up computation...");
    
else
    throw(CORAerror('CORA:noSuitableSolver','SDP'));
end

% ------------------------------ END OF CODE ------------------------------
