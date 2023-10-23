function E = orEllipsoidOA(E)
% orEllipsoidOA - Computes an outer-approximation of the union between
%    ellipsoids
%
% Syntax:
%    E = orEllipsoidOA(E)
%
% Inputs:
%    E - class array of ellipsoid objects
%
% Outputs:
%    E - ellipsoid after union
%
% References:
%   [1] S. Boyd et al. "Convex Optimization"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       15-March-2021
% Last update:   05-July-2022 (VG, remove unecessary input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% collapse cell array
N = length(E);
n = dim(E(1));

% normalize to prevent numerical issues
max_val = max(arrayfun(@(ii)max(svd(E(ii).Q)),1:N));


fac = 0.001;
th = fac*max_val;
th(th==0) = fac;
for i=1:N
    Ei = E(i);
    if ~isFullDim(Ei)
        % if Ei is degenerate, add small pertubation (possible since we
        % compute an overapproximation)
        nd_i = rank(Ei);
        [Ti,Si,~] = svd(Ei.Q);
        si = diag(Si);
        % choose singular value such that rcond not too small
        % bloat to remove degeneracy
        % use max values from other ellipsoids to avoid numerical problems
        % when solving the problem below
        Si = diag([si(1:nd_i);th*ones(n-nd_i,1)]);
        E(i) = ellipsoid(Ti*Si*Ti',Ei.q);
    end
end

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end
persistent isSDPT3
if isempty(isSDPT3)
    isSDPT3 = isSolverInstalled('sdpt3');
end


% find minimum volume ellipsoid spanning union [1]

if isMosek
    % model the problem ourselves using MOSEK
    % MOSEK solves problems of the form 
    %%% (1) min_{x,X}   c'*x + sum_{j=1}^{N} Cj#Xj,
    %%%         s.t.    lb_i <= a^{(i)}'*x+sum_{j=1}^{N}Aij#Xj <= ub_i, i\in{1,...,m},
    %%%                 lb_x <= x <= ub_x,
    %%%                 x \in K,
    %%%                 Xj >=0, j\in {1,...,N}
    % (A#B = sum(sum(A.*B)))
    % The corresponding dual problem is given by
    %%% (2) max_y       -ub'*y_ub + lb'*y_lb - ub_x'*yx_ub + lb_x'*yx_lb,
    %%%         s.t.    c + A'*(y_ub-y_lb) + yx_ub - yx_lb = s,
    %%%                 Cj + sum_{i=1}^{m}Aij*(y_ub-y_lb) >= 0,
    %%%                 s\in dual(K),
    %%%                 y_ub, y_lb, yx_ub, yx_lb: >= 0,
    % where A' = [a^{(1)},a^{(2)},...], and ">=0" means positive
    % semi-definite for matrix variables.
    % Thus, MOSEK needs these values (Cj,Aij,c,A,ub,lb,ub_x,lb_x) and 
    % cone descriptions K
    % For this specific problem, we identify these matrices using the dual
    % problem formulation (easier)

    % in our case, we want to solve the problem
    %%% (3) max_{l,B,bt} det(B)^(1/n),
    %%%         s.t.    [-B+l_j*Qinv_j,     -bt-l_j*q_j,    zeros(n);
    %%%                  (-bt-l_j*q_j)',    1+l_j*c_j,      -bt';
    %%%                  zeros(n),          -bt,            B] >= 0,
    %%%                 l_j >= 0,
    %%%                     j \in {1,...,N},
    %%%   
    % where Qinv_j = inv(E_c{j}.Q), q_j = E_c{j}.q, and c_j =
    % q_j'*Qinv_j*q_j - 1.
    % As this has neither the primal nor the dual form, we reformulate the
    % objective function [3] to eventually get
    %%% (4) max_{...}   z
    %%%         s.t.    [B,Z,Z',diag(diag(Z))] >= 0,
    %%%                 [-B+l_j*Qinv_j,     -bt-l_j*q_j,    zeros(n);
    %%%                  (-bt-l_j*q_j)',    1+l_j*c_j,      -bt';
    %%%                  zeros(n),          -bt,            B] >= 0,
    %%%                 l_j >= 0,
    %%%                     j \in {1,...,N},       
    %%%                 Z = tril(Z),
    %%%                 z <= prod(diag(Z))^(1/n).
    % The last constraint can be modeled using either exponential cones
    % [Sec 5.2, 3], or power cones (see [Sec. 4.2.4, 3]), or rotated 
    % quadratic cones (see [Sec. 3.2.7, 2]).
    % Thus, the problem we need to model is
    %%% (5) max_{...}   z
    %%%         s.t.    [B,Z,Z',diag(diag(Z))] >= 0,
    %%%                 [-B+l_j*Qinv_j,     -bt-l_j*q_j,    zeros(n);
    %%%                  (-bt-l_j*q_j)',    1+l_j*c_j,      -bt';
    %%%                  zeros(n),          -bt,            B] >= 0,
    %%%                 l_j >= 0,
    %%%                     j \in {1,...,N},       
    %%%                 Z = tril(Z),
    %%%                 [diag(Z)',z] in Kb,
    % where Kb is any of the aforementioned cones (for the implementation,
    % this has to be split, and multiple auxiliary variables are
    % introduced).
    % Looking at this problem, it is quiet obvious that modeling this using
    % the dual form is much easier.
    % In this implementation, we use power cones with the splitting as
    % described in [Sec. 4.2.4, 3], where we additionally introduce n-2
    % additional variables to avoid having same variables in different
    % cones (MOSEK does not allow that). We enforce the equality using
    % linear equality constraints.
    
    % first N entries: size of psd matrix for each ellipsoid; last entry:
    % size of psd matrix resulting from reformulation of det(B)^(1/n_nd)
    [~, res] = mosekopt('symbcon echo(0)');
    % Looking at the general dual form, we have N SDP inequalities, and
    % thus need N PSD matrix variables in the primal formulation
    prob.bardim = [(1+2*n)*ones(1,N),2*n];
 
    % initialize Ab
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subk = [];
    prob.bara.subl = [];
    prob.bara.val = [];
    
    gsum = @(x) x.*(x+1)./2;

    % number of variables for n_nd symmetric matrix
    nS = gsum(n);
    % indices for lower triangular part of n_nd symmetric matrix
    [rl_S,cl_S] = find(tril(ones(n))); 
    [ru_S,cu_S] = find(triu(ones(n)));
    indl_S = sub2ind([n,n],rl_S,cl_S);

    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_bt = id_max+(1:n); id_max = max(id_bt);
    id_B = id_max+(1:nS); id_max = max(id_B);
    id_Zt = id_max+(1:nS); id_max = max(id_Zt);
    % ids for diagonal entries (upper triangular matrix)
    ii_diag_u = gsum(1:n);
    id_Zd = id_Zt(ii_diag_u);
    n_t = n-2;

    % same for all i
    prob.barc.subj = 1:N;
    prob.barc.subk = (n+1)*ones(1,N);
    prob.barc.subl = (n+1)*ones(1,N);
    prob.barc.val = ones(1,N);
    
    for i=1:N
        
        % Ab
        % all are negated due to dual sdp form (C_j -sum... A_ij )
        % i: which variable? j: which psd "constraint"? k,l: col,row?
        % values for "l constraints"
        Qinv_i = inv(E(i).Q);
        c_i = E(i).q'*Qinv_i*E(i).q-1;
        a_subi_l_i = id_l(i)*ones(1,nS+n+1);
        a_subj_l_i = i*ones(1,nS+n+1);
        a_subk_l_i = [rl_S',(n+1)*ones(1,n),n+1];
        a_subl_l_i = [cl_S',1:n,n+1];
        a_val_l_i = [-Qinv_i(indl_S)',(Qinv_i*E(i).q)',-c_i];
        
        % values for "bt constraints"
        a_subi_bt_i = [id_bt,id_bt];
        a_subj_bt_i = i*ones(1,2*n);
        a_subk_bt_i = [(n+1)*ones(1,n),n+1+(1:n)];
        a_subl_bt_i = [1:n,(n+1)*ones(1,n)];
        a_val_bt_i = ones(1,2*n);

        % values for "B constraints"
        a_subi_B_i = [id_B,id_B];
        a_subj_B_i = i*ones(1,2*nS);
        a_subk_B_i = [rl_S',n+1+rl_S'];
        a_subl_B_i = [cl_S',n+1+cl_S'];
        a_val_B_i = [ones(1,nS),-ones(1,nS)];

        % collect values
        prob.bara.subi = [prob.bara.subi, a_subi_l_i,a_subi_bt_i,a_subi_B_i];
        prob.bara.subj = [prob.bara.subj, a_subj_l_i,a_subj_bt_i,a_subj_B_i];
        prob.bara.subk = [prob.bara.subk, a_subk_l_i,a_subk_bt_i,a_subk_B_i];
        prob.bara.subl = [prob.bara.subl, a_subl_l_i,a_subl_bt_i,a_subl_B_i];
        prob.bara.val = [prob.bara.val, a_val_l_i,a_val_bt_i,a_val_B_i];
    end
    
    % constraints due to geomean(B) reformulation
    % formulate [B,Z;Z',diag(Z)]>=0
    % Cb
    % nothing to do...

    % Ab
    % entries for B
    a_subi_B_det = id_B;
    a_subj_B_det = (N+1)*ones(1,nS);
    a_subk_B_det = rl_S';
    a_subl_B_det = cl_S';
    a_val_B_det = -ones(1,nS);

    % entries for LT Z & diag(Z)
    a_subi_Z_det = [id_Zt,id_Zd];
    a_subj_Z_det = (N+1)*ones(1,nS+n);
    a_subk_Z_det = [n+ru_S',n+(1:n)];
    a_subl_Z_det = [cu_S',n+(1:n)];
    a_val_Z_det = [-ones(1,nS),-ones(1,n)];

    % collect
    prob.bara.subi = [prob.bara.subi,a_subi_B_det,a_subi_Z_det];
    prob.bara.subj = [prob.bara.subj,a_subj_B_det,a_subj_Z_det];
    prob.bara.subk = [prob.bara.subk,a_subk_B_det,a_subk_Z_det];
    prob.bara.subl = [prob.bara.subl,a_subl_B_det,a_subl_Z_det];
    prob.bara.val = [prob.bara.val,a_val_B_det,a_val_Z_det];
    
    % we need 1 additional variable for z, and then sum_{i=1}{l-1}2^(l-k)
    % additional variables (where l = ceil(log2(n_nd))) (only required if
    % there is a lower level, e.g. )
    l = ceil(log2(n));
    % we for first level, there are 2^(l-1) aux vars needed, then 2^(l-2)
    % etc. Last level need id_z
    n_t = sum(2.^(l-(1:l-1)));
    id_t = max(id_Zt)+(1:n_t);
    if isempty(id_t)
        id_z = max(id_Zt+1);
    else
        id_z = max(id_t)+1;
    end
    id_t_ext = [id_t,id_z];
    % first, fill in upper-most level, until we "run out of" id_Zd
    % fill in with additional idz if not even 
    id_Zd_ext = id_Zd;
    id_y = [id_l,id_bt,id_B,id_Zt,id_t,id_z];
    n_y = length(id_y);
    if mod(n,2)~=0
        % not even
        id_Zd_ext = [id_Zd,id_z];
    end
    % number of cones
    n_cones = sum(2.^(l-(1:l)));
    At = zeros(3*n_cones,n_y);
    for i=1:ceil(n/2)
        At(3*(i-1)+1:3*i,[id_Zd_ext([2*(i-1)+1,2*i]),id_t_ext(i)]) = ...
                                [0.5,0.5,0;0.5,-0.5,0;0,0,1];
    end
    % cones which are (z,0,t_i) (because both z_ii and z_jj are = z)
    for i=ceil(n/2)+1:2^(l-1)
        At(3*(i-1)+1:3*i,[id_z,id_t(i)]) = [1,0;0,0;0,1];
    end
    % first level done: since we have sum(2.^(l-(1:l))) cones, there are
    % ... remaining. Fill in all remaining ones 
    ii_t_prev = 1;
    for i=2^(l-1)+1:n_cones
        At(3*(i-1)+1:3*i,[id_t_ext([ii_t_prev,ii_t_prev+1]),id_t_ext(i)]) = ...
                                [0.5,0.5,0;0.5,-0.5,0;0,0,1];
        ii_t_prev = ii_t_prev + 2;
    end
    % by additionally introducing rows [-eye(N),zeros(N,n_y-N)], we
    % "produce" N additional primal scalar variables. By then setting
    % prob.blx = [-inf(size(At,1),1);zeros(N,1)], prob.bux =
    % inf(size(At,1)+N,1), the last N rows of the new At, i.e.
    %%%-At*y + yx_ub - yx_lb = s ,
    % are
    %%% l + 0 -yx_lb_N  = 0 <=> l = yx_lb_N >= 0,
    % where yx_lb_N are the dual variables of the state bound of the newly
    % introduced N state variables.
    At = [At;eye(N),zeros(N,n_y-N)];
    At = -At;

    % convert to sparse matrix
    prob.a = sparse(At');

    % set state bounds accordingly
    prob.blx = [-inf(1,size(At,1)-N),zeros(1,N)];
    prob.bux = inf(1,size(At,1));

    % objective in scalar variables
    prob.c = zeros(1,size(At,1));

    % model cones using quadratic cones (see Mosek cookbox v2 for details)
    prob.cones.type = res.symbcon.MSK_CT_QUAD*ones(1,n_cones);
    prob.cones.sub = 1:3*n_cones;
    prob.cones.subptr = 1+3*(0:n_cones-1);
    
    % inequality bounds (equality constraints just have same lb and ub)
    prob.blc = [zeros(1,n_y-1),1];
    prob.buc = prob.blc;

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
                E = ellipsoid;
                return;
            else
                % this might be a little too harsh
                throw(CORAerror('CORA:solverIssue','mosek'));
            end
        end
    else
        % might be a little harsh
        throw(CORAerror('CORA:solverIssue','mosek'));
    end

    
    % get dual solution vector
    % since y = y_lb - y_ub
    y_sol = res.sol.itr.slc-res.sol.itr.suc;
    
    % extract the matrix B
    ind_B_vec = ismember(id_y,id_B);
    B_sol_vec = y_sol(ind_B_vec);
    B_sol = zeros(n);
    B_sol(tril(ones(n))==1) = B_sol_vec;
    B_sol = B_sol + tril(B_sol,-1)';
    
    % extract bt vector
    ind_bt = ismember(id_y,id_bt);
    bt_sol = y_sol(ind_bt);

    % construct ellipsoid parameters
    A2 = B_sol;
    q = -A2\bt_sol;
    Q = inv(A2);
   
elseif isSDPT3
    % IMPORTANT: Vectorization is upper-triangular (in contrast to MOSEK,
    % which is lower-triangular)
    % ALSO: For some reason, off-diagonal elements when stacking upper- 
    % triangular matrices should be multiplied by sqrt(2)...
    
    gsum = @(x) x.*(x+1)./2;

    % number of variables for n_nd symmetric matrix
    nS = gsum(n);
    % indices for lower triangular part of n_nd symmetric matrix
    [ru_S,cu_S] = find(triu(ones(n))); 
    indu_S = sub2ind([n,n],ru_S,cu_S);
    iMask_off = ~ismember(1:nS,gsum(1:n));
    
    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_bt = id_max+(1:n); id_max = max(id_bt);
    id_B = id_max+(1:nS);
    id_B_full = zeros(n);
    id_B_full(indu_S) = id_B;
    id_B_full = id_B_full + triu(id_B_full,1)';

    blk = cell(1+N+1,2);
    blk(1,:) = {'l',N};
    blk(2:end,1) = {'s'};
    % psd matrices for each of the N ellipsoids
    blk(1+(1:N),2) = {1+2*n};
    % psd matrix B
    blk(1+N+1,2) = {n};
    

    C = cell(1+N+1,1);
    At = cell(1+N+1,1);
    
    % y = [l;bt;Bu_vec]
    n_y = N+n+nS;

    nM = gsum(1+2*n);
    id_M_vec = 1:nM;
    id_M = zeros(1+2*n);
    id_M(triu(ones(1+2*n))~=0) = id_M_vec;

    % first linear block to ensure l>=0
    C{1} = sparse(N,1);

    At{1} = sparse(1:N,1:N,-ones(1,N),N,n_y);

    % constraint matrices
    % bt
    ind_bt_t = sub2ind(size(id_M),(n+1)*ones(1,n),n+1+(1:n));
    idMu_bt_t = id_M(ind_bt_t);
    % get index for upper triangular matrix
    a_subk_bt = [nS+(1:n),idMu_bt_t];
    a_subl_bt = [id_bt,id_bt];
    a_val_bt = sqrt(2)*ones(1,2*n);

    % B
    % get index for upper triangular matrix
    ind_B_2 = sub2ind(size(id_M),n+1+ru_S',n+1+cu_S');
    idBu_2 = id_M(ind_B_2);
    a_subk_B = [1:nS,idBu_2];
    a_subl_B = [id_B,id_B];
    % scale off-diag elements with sqrt(2)
    ind_B_off = setdiff(1:nS,gsum(1:n));
    a_val_b1 = ones(1,nS);
    a_val_b1(ind_B_off) = sqrt(2)*a_val_b1(ind_B_off);
    a_val_b2 = -ones(1,nS);
    a_val_b2(ind_B_off) = sqrt(2)*a_val_b2(ind_B_off);
    a_val_B = [a_val_b1,a_val_b2];

    for i=1:N
        % objective matrix
        c_subk_i = n+1;
        c_subl_i = n+1;
        c_val_i = 1;
        C_l_i = sparse(c_subk_i,c_subl_i,c_val_i,blk{1+i,2},blk{1+i,2});
        C{1+i} = C_l_i + triu(C_l_i,1)';
        
        % li
        Qi_inv = inv(E(i).Q);
        ci = E(i).q'*Qi_inv*E(i).q-1;
        
        a_subk_li = [1:nS,nS+(1:n),nS+n+1];
        a_subl_li = id_l(i)*ones(1,nS+n+1);
        Qi_inv_vec = Qi_inv(indu_S)';
        Qi_inv_vec(iMask_off) = sqrt(2)*Qi_inv_vec(iMask_off);
        a_val_li = [-Qi_inv_vec,sqrt(2)*(Qi_inv*E(i).q)',-ci];

        % constraint matrices
        At{1+i} = sparse([a_subk_li,a_subk_bt,a_subk_B],...
                    [a_subl_li,a_subl_bt,a_subl_B],...
                    [a_val_li,a_val_bt,a_val_B],nM,n_y);
    end
    
    % objective matrix for B>=0
    C{1+N+1} = sparse(n,n);

    % constraint matrix for B>=0 (ensuring that "psd-B" and "y-B" are
    % equal)
    aend_subk_B = 1:nS;
    aend_subl_B = id_B;
    % extract diagonal and non-diagonal entries to scale non-diagonal
    % entries with sqrt(2)
    ind_diag_B = ismember(id_B,diag(id_B_full)');
    aend_val_B(ind_diag_B) = -ones(1,sum(ind_diag_B));
    aend_val_B(~ind_diag_B) = -sqrt(2)*ones(1,sum(~ind_diag_B));
    At{1+N+1} = sparse(aend_subk_B,aend_subl_B,aend_val_B,nS,n_y);

    % set logdet terms
    OPTIONS.parbarrier = [repmat({0},1+N,1);{1}];
    OPTIONS.printlevel = 0;

    % since sdpt3 supports logdet explicitly, b = 0
    b = sparse(n_y,1);

    % call solver
    [~,~,y,~,info] = sqlp(blk,At,C,b,OPTIONS);
    
    if info.termcode~=0
        throw(CORAerror('CORA:solverIssue','sdpt3'));
    end

    % extract solution B and d
    % y = [l;bt;Bu_vec]
    bt_sol = y(N+(1:n));
    Bu_vec = y(N+n+(1:nS));
    B_sol = zeros(n);
    B_sol(triu(ones(n))==1) = Bu_vec;
    B_sol = B_sol + triu(B_sol,1)';

    % construct ellipsoid parameters
    A2 = B_sol;
    q = -A2\bt_sol;
    Q = inv(A2);

elseif isYalmipInstalled()
    A2 = sdpvar(n); % pd-ness of A ensured by yalmip internally
    bt = sdpvar(n,1);
    l = sdpvar(N,1);
    f_obj = -1/2*geomean(A2);
    C = [];
    for i=1:N
        Qiinv = inv(E(i).Q);
        qi = E(i).q;
        Ci = [A2-l(i)*Qiinv,    bt+l(i)*Qiinv*qi,            zeros(n);
              (bt+l(i)*Qiinv*qi)',-1-l(i)*(qi'*Qiinv*qi-1),  bt';
              zeros(n),          bt,                         -A2];
        C = [C, Ci<=0];
    end
    opts = sdpsettings;
    opts.verbose = 0;
    sol = optimize([C,l>=0],f_obj,opts);
    warning("YALMIP was used to model the problem - " + ...
        "consider installing a supported solver to speed up computation...");
    if any(sol.problem == [-2,-3,-4])
        throw(CORAerror('CORA:YALMIP',...
            'This function requires an SDP solver compatible with YALMIP!'));
    elseif sol.problem ~= 0
        throw(CORAerror('CORA:solverIssue'));
    end
    
    % extract ellipsoid
    Q = inv(value(A2));
    A = sqrtm(value(A2));
    b = A\value(bt);
    q = -Q*A*b;
    
else
    throw(CORAerror('CORA:noSuitableSolver','SDP'));
end

% backtransform
E = ellipsoid(Q,q);

% ------------------------------ END OF CODE ------------------------------
