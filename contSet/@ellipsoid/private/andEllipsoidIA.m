function E = andEllipsoidIA(E,E2)
% andEllipsoidIA - Computes the inner-approximation of the intersection
%    between multiple ellipsoids
%
% Syntax:
%    E = andEllipsoidIA(E,E2)
%
% Inputs:
%    E - ellipsoid object
%    E_c - ellipsoid object or cell array of ellipsoid objects
%
% Outputs:
%    E - ellipsoid after intersection
%
% References:
%   [1] S. Boyd et al. "Convex Optimization"
%   [2] MOSEK Modelling Cookbook v2
%   (https://docs.mosek.com/MOSEKModelingCookbook-v2.pdf)
%   [3] MOSEK Modelling Cookbook v3
%   (https://docs.mosek.com/MOSEKModelingCookbook-a4paper.pdf)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: andEllipsoidOA

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   13-May-2022
%                04-July-2022 (VG, convert cell array to class array)
%                06-September-2022 (VG, resolved bug for 1D-case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

E = [E;E2(:)];
ind_dMask = ~isFullDim(E);
if sum(ind_dMask)>1
    throw(CORAerror('CORA:degenerateSet',...
        'At most one ellipsoid can be degenerate!'));
end

% make at most last array element degenerate
Ed = E(ind_dMask);
End = E(~ind_dMask);
E = [End;Ed];

n = dim(E(1));
N = length(E);
I = eye(n);
T = eye(n);
x_rem = zeros(0,1);

if ~isFullDim(E(end))
    [T,~,~] = svd(E(end).Q);
    nt = rank(E(end));
    % transform all ellipsoids
    E = T'*E;
    % project
    x_rem = E(end).q(nt+1:end);
    E(end) = project(E(end),1:nt);
    % resolve x_rem in E(1:end-1) by cutting with hyperplanes 
    % I(nt+1:end,:)*xt = x_rem
    for i=1:N-1
        % intersect each cell element with all hyperplanes
        for j=1:(n-nt)
            Hi = conHyperplane(I(nt+j,:),x_rem(j));
            E(i) = and_(E(i),Hi,'outer');
            % check if intersection is empty; if yes, overall result is
            % empty
            if representsa_(E(i),'emptySet',eps)
                E = ellipsoid;
                return;
            end
        end
        % E_c{i} now also has zeros at nt+1:n
        E(i) = project(E(i),1:nt);
    end
end

n_nd = E(1).dim;

% if 1d remaining, use interval arithmetic
if n_nd==1
    Ires = interval(E(1));
    for j=2:N
        Ires = and_(Ires,interval(E(j)),'exact');
    end
    if representsa_(Ires,'emptySet',eps)
        E = ellipsoid;
        return;
    end
    qt = center(Ires);
    Qt = rad(Ires)^2;

    % backtransform
    E = T*ellipsoid([Qt,zeros(n_nd,n-n_nd);zeros(n-n_nd,n)],[qt;x_rem]);
    return
end

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
    % MOSEK solves problems of the form (var bounds omitted, not needed for
    % us)
    %%% (1) min_{x,X}   c'*x + sum_{j=1}^{N} Cj#Xj,
    %%%         s.t.    lb_i <= a^{(i)}'*x+sum_{j=1}{N}Aij#Xj <= ub_i, i\in{1,...,m},
    %%%                 x \in K,
    %%%                 Xj >=0, j\in {1,...,N}
    % (A#B = sum(sum(A.*B)))
    % The corresponding dual problem is given by
    %%% (2) max_y       -ub'*y_ub + lb'*y_lb - ub_x'*yx_ub + lb_x'*yx_lb,
    %%%         s.t.    c + A'*(y_ub-y_lb) yx_ub - yx_lb = s,
    %%%                 Cj + sum_{i=1}{m}Aij*(y_ub-y_lb) >=0
    %%%                 s\in dual(K),
    %%%                 y_ub>=0, y_lb >=0, yx_ub>=0, yx_lb>=0,
    % where A' = [a^{(1)},a^{(2)},...], and ">=0" means positive
    % semi-definite for matrix variables.
    % Thus, MOSEK needs these matrices (Cj, Aij, c, A) and cone
    % descriptions K
    % For this specific problem, we identify these matrices using the dual
    % problem formulation (easier)

    % in our case, we want to solve the problem
    %%% (3) max_{B,d,l} det(B)^(1/n_nd),
    %%%         s.t.    [-l_j+1,        zeros(1,n_nd),  (d-q^{(j)})';
    %%%                  zeros(n_nd,1), l_j*eye(n_nd),  B;
    %%%                  d-q^{(j)},     B,              Q^{(j)}]>=0,
    %%%                     j \in {1,...,N},
    %%%                 l_j >= 0
    % As this has neither the primal nor the dual form, we reformulate the
    % objective function [3] to eventually get
    %%% (4) max_{...}   z
    %%%         s.t.    [B,Z,Z',diag(diag(Z))] >= 0,
    %%%                 [-l_j+1,        zeros(1,n_nd),  (d-q^{(j)})';
    %%%                  zeros(n_nd,1), l_j*eye(n_nd),  B;
    %%%                  d-q^{(j)},     B,              Q^{(j)}]>=0,
    %%%                     j \in {1,...,N},       
    %%%                 Z = tril(Z),
    %%%                 z <= prod(diag(Z))^(1/n_nd),
    %%%                 l_j >= 0
    % The last constraint can be modeled using either exponential cones
    % [Sec 5.2, 3], or power cones (see [Sec. 4.2.4, 3]), or rotated 
    % quadratic cones (see [Sec. 3.2.7, 2]).
    % Thus, the problem we need to model is
    %%% (5) max_{...}   z
    %%%         s.t.    [B,Z,Z',diag(diag(Z))] >= 0,
    %%%                 [-l_j+1,        zeros(1,n_nd),  (d-q^{(j)})';
    %%%                  zeros(n_nd,1), l_j*eye(n_nd),  B;
    %%%                  d-q^{(j)},     B,              Q^{(j)}]>=0,
    %%%                     j \in {1,...,N},       
    %%%                 Z = tril(Z),
    %%%                 [diag(Z)',z] in Kb,
    %%%                 l_j>=0,
    % where Kb is any of the aforementioned cones (for the implementation,
    % this has to be split, and multiple auxiliary variables are
    % introduced).
    % Looking at this problem, it is quiet obvious that modeling this using
    % the dual form is much easier.
    % In this implementation, we use power cones with the splitting as
    % described in [Sec. 4.2.4, 3], where we additionally introduce n_nd-2
    % additional variables to avoid having same variables in different
    % cones (MOSEK does not allow that). We enforce the equality using
    % linear equality constraints.
    % Lastly, here we only require a simplified version of the dual problem
    % with lb = ub (equality constraint), which means we can introduce a
    % new variable y = y_lb - y_ub (in \mathbb{R}^m), which transforms (2)
    % into (with b = lb = ub)
    %%% (6) max_y       b'*y,
    %%%         s.t.    c - A'*y = s,
    %%%                 Cj - sum_{i=1}{m}Aij*y >=0
    %%%                 s\in dual(K).
    

    % first N entries: size of psd matrix for each ellipsoid; last entry:
    % size of psd matrix resulting from reformulation of det(B)^(1/n_nd)
    [~, res] = mosekopt('symbcon echo(0)');
    % Looking at the general dual form, we have N SDP inequalities, and
    % thus need N PSD matrix variables in the primal formulation
    prob.bardim = [(1+2*n_nd)*ones(1,N),2*n_nd];
    
    % initialize Cb
    prob.barc.subj = [];
    prob.barc.subk = [];
    prob.barc.subl = [];
    prob.barc.val = [];
    % initialize Ab
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subk = [];
    prob.bara.subl = [];
    prob.bara.val = [];
    
    gsum = @(x) x.*(x+1)./2;

    % number of variables for n_nd symmetric matrix
    nS = gsum(n_nd);
    % indices for lower triangular part of n_nd symmetric matrix
    [rl_S,cl_S] = find(tril(ones(n_nd))); 
    [ru_S,cu_S] = find(triu(ones(n_nd)));
    indl_S = sub2ind([n_nd,n_nd],rl_S,cl_S);

    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_d = id_max+(1:n_nd); id_max = max(id_d);
    id_B = id_max+(1:nS); id_max = max(id_B);
    id_Zt = id_max+(1:nS); id_max = max(id_Zt);
    % ids for diagonal entries (upper triangular matrix)
    ii_diag_u = gsum(1:n_nd);
    id_Zd = id_Zt(ii_diag_u);

    % this gives us
    id_B_full = zeros(n_nd);
    id_B_full(indl_S) = id_B;
    id_B_full = id_B_full + tril(id_B_full,-1)';
    id_B_vec_full = id_B_full(:)';
    [r_B_full,c_B_full] = ind2sub([n_nd,n_nd],1:n_nd^2);
    
    for i=1:N
        % Cb = [1,0',-E_c{i}.q';0;-E_c{i}.q,0,E_c{i}.Q]
        c_subj_i = i*ones(1,1+n_nd+nS);
        c_subk_i = [1,1+n_nd+(1:n_nd),1+n_nd+rl_S'];
        c_subl_i = [1,ones(1,n_nd),1+n_nd+cl_S'];
        c_val_i = [1,-E(i).q',E(i).Q(indl_S)'];
        
        prob.barc.subj = [prob.barc.subj,c_subj_i];
        prob.barc.subk = [prob.barc.subk,c_subk_i];
        prob.barc.subl = [prob.barc.subl,c_subl_i];
        prob.barc.val = [prob.barc.val,c_val_i];

        % Ab
        % all are negated due to dual sdp form (C_j -sum... A_ij )
        
        % values for "l constraints"
        a_subi_l_i = id_l(i)*ones(1,1+n_nd);
        a_subj_l_i = i*ones(1,1+n_nd);
        a_subk_l_i = [1,1+(1:n_nd)];
        a_subl_l_i = [1,1+(1:n_nd)];
        a_val_l_i = [1,-ones(1,n_nd)];
        
        % values for "d constraints"
        a_subi_d_i = id_d;
        a_subj_d_i = i*ones(1,n_nd);
        a_subk_d_i = 1+n_nd+(1:n_nd);
        a_subl_d_i = ones(1,n_nd);
        a_val_d_i = -ones(1,n_nd);

        % values for "B constraints"
        a_subi_B_i = id_B_vec_full;
        a_subj_B_i = i*ones(1,n_nd^2);
        a_subk_B_i = 1+n_nd+r_B_full;
        a_subl_B_i = 1+c_B_full;
        a_val_B_i = -ones(1,n_nd^2);

        % collect values
        prob.bara.subi = [prob.bara.subi, a_subi_l_i,a_subi_d_i,a_subi_B_i];
        prob.bara.subj = [prob.bara.subj, a_subj_l_i,a_subj_d_i,a_subj_B_i];
        prob.bara.subk = [prob.bara.subk, a_subk_l_i,a_subk_d_i,a_subk_B_i];
        prob.bara.subl = [prob.bara.subl, a_subl_l_i,a_subl_d_i,a_subl_B_i];
        prob.bara.val = [prob.bara.val, a_val_l_i,a_val_d_i,a_val_B_i];
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
    a_subj_Z_det = (N+1)*ones(1,nS+n_nd);
    a_subk_Z_det = [n_nd+ru_S',n_nd+(1:n_nd)];
    a_subl_Z_det = [cu_S',n_nd+(1:n_nd)];
    a_val_Z_det = [-ones(1,nS),-ones(1,n_nd)];

    % collect
    prob.bara.subi = [prob.bara.subi,a_subi_B_det,a_subi_Z_det];
    prob.bara.subj = [prob.bara.subj,a_subj_B_det,a_subj_Z_det];
    prob.bara.subk = [prob.bara.subk,a_subk_B_det,a_subk_Z_det];
    prob.bara.subl = [prob.bara.subl,a_subl_B_det,a_subl_Z_det];
    prob.bara.val = [prob.bara.val,a_val_B_det,a_val_Z_det];
    
    % we need 1 additional variable for z, and then sum_{i=1}{l-1}2^(l-k)
    % additional variables (where l = ceil(log2(n_nd))) (only required if
    % there is a lower level, e.g. )
    l = ceil(log2(n_nd));
    % for first level, there are 2^(l-1) aux vars needed, then 2^(l-2)
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
    id_y = [id_l,id_d,id_B,id_Zt,id_t,id_z];
    n_y = length(id_y);
    if mod(n_nd,2)~=0
        % not even
        id_Zd_ext = [id_Zd,id_z];
    end
    % number of cones
    n_cones = sum(2.^(l-(1:l)));
    At = zeros(3*n_cones,n_y);
    for i=1:ceil(n_nd/2)
        At(3*(i-1)+1:3*i,[id_Zd_ext([2*(i-1)+1,2*i]),id_t_ext(i)]) = ...
                                [0.5,0.5,0;0,0,1;0.5,-0.5,0];
    end
    % cones which are (z,0,t_i) (because both z_ii and z_jj are = z)
    for i=ceil(n_nd/2)+1:2^(l-1)
        At(3*(i-1)+1:3*i,[id_z,id_t(i)]) = [1,0;0,0;0,1];
    end
    % first level done: since we have sum(2.^(l-(1:l))) cones, there are
    % ... remaining. Fill in all remaining ones 
    ii_t_prev = 1;
    for i=2^(l-1)+1:n_cones
        At(3*(i-1)+1:3*i,[id_t_ext([ii_t_prev,ii_t_prev+1]),id_t_ext(i)]) = ...
                                [0.5,0.5,0;0,0,1;0.5,-0.5,0];
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
    At = [eye(N),zeros(N,n_y-N);At];
    At = -At;
    prob.a = sparse(At');

    % set state bounds accordingly
    prob.blx = [zeros(1,N),-inf(1,size(At,1)-N)];
    prob.bux = inf(1,size(At,1));

    % objective in scalar variables
    prob.c = zeros(1,size(At,1));

    % model cones using quadratic cones (see Mosek cookbox v2 for details)
    prob.cones.type = res.symbcon.MSK_CT_QUAD*ones(1,n_cones);
    prob.cones.sub = 2+(1:3*n_cones);
    prob.cones.subptr = 1+3*(0:n_cones-1);

    % inequality constraints (equality constraints just have same lb and ub)
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
        throw(CORAerror('CORA:solverIssue','mosek'));
    end

    
    % get dual solution vector
    % since y = y_lb - y_ub
    y_sol = res.sol.itr.slc-res.sol.itr.suc;
    
    % extract the "shape" matrix B
    ind_B_vec = ismember(id_y,id_B);
    B_sol_vec = y_sol(ind_B_vec);
    B_sol = zeros(n_nd);
    B_sol(tril(ones(n_nd))==1) = B_sol_vec;
    B_sol = B_sol + tril(B_sol,-1)';
    
    % extract d vector
    ind_d = ismember(id_y,id_d);
    d_sol = y_sol(ind_d);

    % construct Qt matrix
    Qt = B_sol^2;
    qt = d_sol;

elseif isSDPT3
    % IMPORTANT: Vectorization is upper-triangular (in contrast to MOSEK,
    % which is lower-triangular)
    % ALSO: For some reason, off-diagonal elements when stacking upper- 
    % triangular matrices should be multiplied by sqrt(2)...
    
    gsum = @(x) x.*(x+1)./2;

    % number of variables for n_nd symmetric matrix
    nS = gsum(n_nd);
    % indices for lower triangular part of n_nd symmetric matrix
    [ru_S,cu_S] = find(triu(ones(n_nd))); 
    indu_S = sub2ind([n_nd,n_nd],ru_S,cu_S);
    
    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_d = id_max+(1:n_nd); id_max = max(id_d);
    id_B = id_max+(1:nS);
    id_B_full = zeros(n_nd);
    id_B_full(indu_S) = id_B;
    id_B_full = id_B_full + triu(id_B_full,1)';
    id_B_vec_full = id_B_full(:)';

    blk = cell(1+N+1,2);
    blk(1,:) = {'l',N};
    blk(2:end,1) = {'s'};
    % psd matrices for each of the N ellipsoids
    blk(1+(1:N),2) = {1+2*n_nd};
    % psd matrix B
    blk(1+N+1,2) = {n_nd};
    
    C = cell(1+N+1,1);
    At = cell(1+N+1,1);
    
    % y = [l;d;Bu_vec]
    n_y = N+n_nd+nS;

    nM = gsum(1+2*n_nd);
    id_M_vec = 1:nM;
    id_M = zeros(1+2*n_nd);
    id_M(triu(ones(1+2*n_nd))~=0) = id_M_vec;
    
    % first linear block to ensure l>=0
    C{1} = sparse(N,1);

    At{1} = sparse(1:N,1:N,-ones(1,N),N,n_y);

    % constraint matrices
    % li
    % get index for upper triangular matrix
    tmp = diag(id_M);
    a_subk_li = tmp(1:1+n_nd)';
    a_Subl_l = (id_l').*ones(N,1+n_nd);
    a_val_li = [1,-ones(1,n_nd)];

    % d
    % get index for upper triangular matrix
    a_subk_d = id_M(1,1+n_nd+(1:n_nd));
    a_subl_d = id_d;
    a_val_d = -sqrt(2)*ones(1,n_nd);

    % B
    % get index for upper triangular matrix
    a_subk_B = reshape(id_M(1+(1:n_nd),1+n_nd+(1:n_nd)),1,[]);
    a_subl_B = id_B_vec_full;
    a_val_B = -sqrt(2)*ones(1,n_nd^2);

    for i=1:N
        % objective matrix
        c_subk_i = [1,ones(1,n_nd),1+n_nd+ru_S'];
        c_subl_i = [1,1+n_nd+(1:n_nd),1+n_nd+cu_S'];
        c_val_i = [1,-E(i).q',E(i).Q(indu_S)'];
        C_l_i = sparse(c_subk_i,c_subl_i,c_val_i,blk{1+i,2},blk{1+i,2});
        C{1+i} = C_l_i + triu(C_l_i,1)';

        % constraint matrices
        At{1+i} = sparse([a_subk_li,a_subk_d,a_subk_B],...
                    [a_Subl_l(i,:),a_subl_d,a_subl_B],...
                    [a_val_li,a_val_d,a_val_B],nM,n_y);
    end
    
    % objective matrix for B>=0
    C{1+N+1} = sparse(n_nd,n_nd);

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
    [~,~,y,Z,info] = sqlp(blk,At,C,b,OPTIONS);
    
    if info.termcode~=0
        % either primal and/or dual infeasible (should only ever occur
        % both), or numerical issues
        if any(info.termcode==[1,2])
            % primal and/or dual infeasible
            E = ellipsoid;
            return;
        end
        % do not know error code, assume some problem occured
        throw(CORAerror('CORA:solverIssue','sdpt3'));
    end

    % extract solution B and d
    % y = [l;d;Bu_vec]
    % could also be extracted from y
    B_sol = Z{1+N+1};
    d_sol = y(N+1:N+n_nd);

    % construct 
    Qt = B_sol^2;
    qt = d_sol;

elseif isYalmipInstalled()
    % formalize and solve optimization problem
    B_sol = sdpvar(n_nd);
    l = sdpvar(N,1);
    d_sol = sdpvar(n_nd,1);
    f_obj = -geomean(B_sol);
    C = [];
    for j=1:N
        Aiinv = E(j).Q;
        bi = -inv(E(j).Q)*E(j).q;
        % formulate constraint
        Ci = [-l(j)+1,zeros(1,n_nd),       (d_sol-E(j).q)';
              zeros(n_nd,1),           l(j)*eye(n_nd),    B_sol;
              d_sol-E(j).q,           B_sol,              Aiinv];
        C = [C,Ci>=0];
    end
    sol = optimize([C,l>=0],f_obj,sdpsettings('verbose',0));
    warning("YALMIP was used to model the problem - " + ...
        "consider installing a supported solver to speed up computation...");
    if sol.problem==0
        Qt = value(B_sol)^2;
        qt = value(d_sol);
    elseif sol.problem==1
        E = ellipsoid;
        return;
    elseif any(sol.problem==[-2,-3,-9])
        throw(CORAerror('CORA:noSuitableSolver','SDP'));
    else
        throw(CORAerror('CORA:solverIssue'));
    end

else
    throw(CORAerror('CORA:noSuitableSolver','SDP'));
end

% backtransform
E = T*ellipsoid([Qt,zeros(n_nd,n-n_nd);zeros(n-n_nd,n)],[qt;x_rem]);

% ------------------------------ END OF CODE ------------------------------
