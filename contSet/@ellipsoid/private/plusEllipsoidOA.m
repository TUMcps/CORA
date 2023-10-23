function E = plusEllipsoidOA(E)
% plusEllipsoidOA - Computes the smallest-volume outer-approximation of the
%    Minkowski sum of ellipsoids
%
% Syntax:
%    E = plusEllipsoidOA(E)
%
% Inputs:
%    E - array of ellipsoid objects
%
% Outputs:
%    E - ellipsoid object
%
% References:
%   [1] Boyd, Stephen, et al. "Linear matrix inequalities in system and 
%       control theory"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   05-July-2022 (VG, class array support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no error checking, already done in parent function

% number of ellipsoids in cell array
N = length(E);

% addin' nothin' to somethin' is still that somethin'
if N == 1
    return
end

% dimension
n = dim(E(1));

% it can happen that we have degenerate ellipsoids (possibly even all;
% their sum does form a non-degenerate ellipsoid, though). Thus, since we
% compute inverses later on, we simply bloat the ellipsoids here.
for i=1:N
    if ~isFullDim(E(i))
        % bloat ellipsoid in appropriate dimensions
        [U,S,~] = svd(E(i).Q);
        nt = rank(E(i));
        s = diag(S);
        s(nt+1:end) = 2*max(s(1)*E(i).TOL,E(i).TOL);
        Q_i = U*diag(s)*U';
        E(i).Q = 1/2*(Q_i+Q_i)';
    end
end

E_c = cell(N,1);
E0 = zeros(n,n*N);
At_c = cell(N,1);
bt_c = cell(N,1);
c_c = cell(N,1);
for i=1:N
    % see [1] for details
    E_c{i} = [zeros(n,(i-1)*n),eye(n),zeros(n,(N-i)*n)];
    E0 = E0 + E_c{i};
    b_i = -E(i).Q\E(i).q;
    At_c{i} = E_c{i}'*(E(i).Q\E_c{i});
    bt_c{i} = E_c{i}'*b_i;
    c_c{i} = E(i).q'*(E(i).Q\E(i).q)-1;
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
    % For details, see "orEllipsoidOA.m" (very similiar)
    
    [~, res] = mosekopt('symbcon echo(0)');

    prob.bardim = [n*N+1+n,2*n];
 
    % initialize Ab
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subk = [];
    prob.bara.subl = [];
    prob.bara.val = [];
    
    gsum = @(x) x.*(x+1)./2;

    % number of variables for n symmetric matrix
    nS = gsum(n);
    % indices for lower triangular part of n symmetric matrix
    [rl_S,cl_S] = find(tril(ones(n))); 
    [ru_S,cu_S] = find(triu(ones(n)));
    indl_S = sub2ind([n,n],rl_S,cl_S);

    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_b = id_max+(1:n); id_max = max(id_b);
    id_B = id_max+(1:nS); id_max = max(id_B);
    id_Zt = id_max+(1:nS); id_max = max(id_Zt);
    % ids for diagonal entries (upper triangular matrix)
    ii_diag_u = gsum(1:n);
    id_Zd = id_Zt(ii_diag_u);

    id_B_full = zeros(n);
    id_B_full(indl_S) = id_B;
    id_B_full = id_B_full + tril(id_B_full,-1)';

    % same for all i
    prob.barc.subj = 1;
    prob.barc.subk = n*N+1;
    prob.barc.subl = n*N+1;
    prob.barc.val = 1;

    % constraint matrix
    % B
    [rl_BB,cl_BB,val_BB] = find(tril(E0'*id_B_full*E0));
    a_subi_B = [val_BB',id_B];
    a_subj_B = ones(1,length(val_BB)+nS);
    a_subk_B = [rl_BB',n*N+1+rl_S'];
    a_subl_B = [cl_BB',n*N+1+cl_S'];
    a_val_B = [ones(1,length(val_BB)),-ones(1,nS)];

    % b
    [r_bb,c_bb,val_bb] = find(id_b*E0);
    a_subi_b = [val_bb,id_b];
    a_subj_b = ones(1,length(val_bb)+n);
    a_subk_b = [n*N+r_bb,n*N+1+(1:n)];
    a_subl_b = [c_bb,(n*N+1)*ones(1,n)];
    a_val_b = ones(1,length(val_bb)+n);
    
    % already collect these values
    prob.bara.subi = [a_subi_B,a_subi_b];
    prob.bara.subj = [a_subj_B,a_subj_b];
    prob.bara.subk = [a_subk_B,a_subk_b];
    prob.bara.subl = [a_subl_B,a_subl_b];
    prob.bara.val = [a_val_B,a_val_b];

    for i=1:N
        [rl_Ati,cl_Ati,val_Ati] = find(tril(At_c{i}));
        n_Ati = length(val_Ati);
        [r_bti_t,c_bti_t,val_bti_t] = find(bt_c{i}');
        n_bti_t = length(val_bti_t);
        % l_i
        a_subi_li = id_l(i)*ones(1,n_Ati+n_bti_t+1);
        a_subj_li = ones(1,n_Ati+n_bti_t+1);
        a_subk_li = [rl_Ati',n*N+r_bti_t,n*N+1];
        a_subl_li = [cl_Ati',c_bti_t,n*N+1];
        a_val_li = [-val_Ati',-val_bti_t,-c_c{i}];

        % collect values
        prob.bara.subi = [prob.bara.subi, a_subi_li];
        prob.bara.subj = [prob.bara.subj, a_subj_li];
        prob.bara.subk = [prob.bara.subk, a_subk_li];
        prob.bara.subl = [prob.bara.subl, a_subl_li];
        prob.bara.val = [prob.bara.val, a_val_li];
    end
    
    % constraints due to geomean(B) reformulation
    % formulate [B,Z;Z',diag(Z)]>=0
    % Cb
    % nothing to do...

    % Ab
    % entries for B
    a_subi_B_det = id_B;
    a_subj_B_det = 2*ones(1,nS);
    a_subk_B_det = rl_S';
    a_subl_B_det = cl_S';
    a_val_B_det = -ones(1,nS);

    % entries for LT Z & diag(Z)
    a_subi_Z_det = [id_Zt,id_Zd];
    a_subj_Z_det = 2*ones(1,nS+n);
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
    id_y = [id_l,id_b,id_B,id_Zt,id_t,id_z];
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
    if res.rcode~=0 || ~strcmp(res.sol.itr.solsta,'OPTIMAL')
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
    ind_b = ismember(id_y,id_b);
    b_sol = y_sol(ind_b);

    % construct ellipsoid parameters
    Q = inv(B_sol);
    q = -Q*value(b_sol);
    E = ellipsoid(Q,q);
   
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
    
    % generate ids to identify same constraints: For our problem
    id_max = 0;
    id_l = id_max+(1:N); id_max = max(id_l);
    id_b = id_max+(1:n); id_max = max(id_b);
    id_B = id_max+(1:nS);
    id_B_full = zeros(n);
    id_B_full(indu_S) = id_B;
    id_B_full = id_B_full + triu(id_B_full,1)';
    
    blk = cell(3,2);
    % one linear block to ensure l>=0, and one psd block for matrix
    % inequality, and one psd block for logdet term
    blk(1,:) = {'l',N};
    blk(2,:) = {'s',n*N+1+n};
    blk(3,:) = {'s',n};

    C = cell(3,1);
    At = cell(3,1);
    
    % y = [l;b;Bu_vec]
    n_y = N+n+nS;

    id_S_vec = 1:nS;
    id_S = zeros(n);
    id_S(triu(ones(n))~=0) = id_S_vec;

    % first linear block to ensure l>=0
    C{1} = sparse(N,1);

    % At_l*y + z_l = c_l
    At{1} = sparse(1:N,1:N,-ones(1,N),N,n_y);

    % constraint matrices
    % B
    [~,~,val_BB] = find(triu(E0'*id_B_full*E0));
    res_c = svec({'s',n*N},{E0'*ones(n)*E0});
    [a_subk_BBtop_,~,a_val_BBtop_] = find(res_c{1});
    % get indices of lower-right part of matrix (-B part)
    % first row of -B
    r_Bl = n*N+1+1;
    nn = n*N+1+n;
    % construct indices of matrix without explicitely constructing the
    % large n*N+1+n matrix
    % only upper triangular indices are accurate
    id_MB = id_S + (r_Bl-1)*(0:nn-r_Bl) + gsum(r_Bl)-1;
    id_MB_u = id_MB(indu_S)';
    a_subk_B = [a_subk_BBtop_',id_MB_u];
    a_subl_B = [val_BB',id_B];
    % split bot into diag and non-diag part
    a_val_Bbot = ones(1,nS);
    % scale off-diag entries by sqrt(2)
    tmp = gsum(1:n*N+1+n);
    iMask_B_od = ~ismember(id_MB_u,tmp(n*N+1+1:end));
    a_val_Bbot(iMask_B_od) = sqrt(2)*a_val_Bbot(iMask_B_od);
    a_val_B = [a_val_BBtop_',-a_val_Bbot];

    % b
    [~,~,val_bb] = find(E0'*id_b');
    % indices of b'
    a_subl_bt = gsum(n*N+1)+gsum(n*N+1:nn-1)-gsum(n*N+1-1);
    a_subk_b = [gsum(n*N)+(1:n*N),a_subl_bt];
    a_subl_b = [val_bb',id_b];
    a_val_b = sqrt(2)*ones(1,n*N+n);
    
    % objective matrix
    % same for all i
    C{2} = sparse(n*N+1,n*N+1,1,n*N+1+n,n*N+1+n);
    
    a_subk = [a_subk_b,a_subk_B];
    a_subl = [a_subl_b,a_subl_B];
    a_val = [a_val_b,a_val_B];

    for i=1:N
        
        [rl_Ati,cl_Ati,val_Ati] = find(tril(At_c{i}));
        [r_bti_t,c_bti_t,val_bti_t] = find(bt_c{i}');
        % l_i

        a_subk_li_ = [rl_Ati',n*N+r_bti_t,n*N+1];
        a_subl_li_ = [cl_Ati',c_bti_t,n*N+1];
        a_val_li_ = [-val_Ati',-val_bti_t,-c_c{i}];
        At_i = sparse(a_subk_li_,a_subl_li_,a_val_li_,n*N+1+n,n*N+1+n)';
        At_i_c = svec(blk(2,:),{At_i});

        [at_subk_li,~,at_val_li] = find(At_i_c{1});
        a_subk_li = at_subk_li';
        a_subl_li = id_l(i)*ones(1,length(a_subk_li));
        a_val_li = at_val_li';
        

        % constraint matrices
        a_subk = [a_subk,a_subk_li];
        a_subl = [a_subl,a_subl_li];
        a_val = [a_val,a_val_li];
    end
    
    At{2} = sparse(a_subk,a_subl,a_val,gsum(n*N+1+n),n_y);

    % objective matrix for B>=0
    C{3} = sparse(n,n);

    % constraint matrix for B>=0 (ensuring that "psd-B" and "y-B" are
    % equal)
    aend_subk_B = 1:nS;
    aend_subl_B = id_B;
    % extract diagonal and non-diagonal entries to scale non-diagonal
    % entries with sqrt(2)
    ind_diag_B = ismember(id_B,diag(id_B_full)');
    aend_val_B(ind_diag_B) = -ones(1,sum(ind_diag_B));
    aend_val_B(~ind_diag_B) = -sqrt(2)*ones(1,sum(~ind_diag_B));
    At{3} = sparse(aend_subk_B,aend_subl_B,aend_val_B,nS,n_y);

    % set logdet terms
    OPTIONS.parbarrier = {0,0,1};
    OPTIONS.printlevel = 0;

    % since sdpt3 supports logdet explicitly, b = 0
    b = sparse(n_y,1);

    % call solver
    [~,~,y,~,info] = sqlp(blk,At,C,b,OPTIONS);
    
    if info.termcode~=0
        throw(CORAerror('CORA:solverIssue','sdpt3'));
    end

    % extract solution B and d
    % y = [l;b;Bu_vec]
    b_sol = y(N+(1:n));
    Bu_vec = y(N+n+(1:nS));
    B_sol = zeros(n);
    B_sol(triu(ones(n))==1) = Bu_vec;
    B_sol = B_sol + triu(B_sol,1)';
    
    % construct ellipsoid parameters
    Q = inv(B_sol);
    q = -Q*value(b_sol);
    E = ellipsoid(Q,q);

elseif isYalmipInstalled()
    B = sdpvar(n);
    b = sdpvar(n,1);
    l = sdpvar(N,1);

    C = [E0'*B*E0,E0'*b,zeros(n*N,n);
         b'*E0,-1,b';
         zeros(n,n*N),b,-B];
    for i=1:N
        C = C - l(i)*[At_c{i},bt_c{i},zeros(n*N,n);
                   bt_c{i}',c_c{i},zeros(1,n);
                   zeros(n,n*N+1+n)];
    end
    sol = optimize([l>=0,C<=0],-logdet(B),sdpsettings('verbose',0));
    warning("YALMIP was used to model the problem - " + ...
        "consider installing a supported solver to speed up computation...");

    % instantiate resulting ellipsoid object
    if sol.problem==0
        Q = inv(value(B));
        q = -Q*value(b);
        E = ellipsoid(Q,q);
    elseif any(sol.problem==[-2,-3,-9])
        throw(CORAerror('CORA:noSuitableSolver','SDP'));
    else
        throw(CORAerror('CORA:solverIssue'));
    end
    
else
    throw(CORAerror('CORA:noSuitableSolver','SDP'));
end

% ------------------------------ END OF CODE ------------------------------
