function E = ellipsoid(obj,varargin)
% ellipsoid - Overapproximates a mptPolytope by an ellipsoid
%
% Syntax:  
%    E = ellipsoid(obj,mode)
%
% Inputs:
%    Z - ellipsoid object
%    mode - (optional): 
%               :'inner' (inner approx)
%               :'outer' (outer approx)
%               :'outer:min-vol' (min-vol outer approx)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      15-March-2021
% Last update:  05-July-2022 (VG: included direct solver support for SDP)
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
mode = setDefaultValues({'outer'},varargin);

% check input arguments
inputArgsCheck({{obj,'att','mptPolytope','scalar'};
                {mode,'str',{'outer','outer:min-vol','inner'}}});


V = vertices(obj);
% remove bias
m = mean(V,2);
V = V-m;

if contains(mode,'outer')
    mm = 'cov';
    if strcmp(mode,'outer:min-vol')
        mm = 'min-vol';
    end
    E = ellipsoid.enclosePoints(V,mm);
else
    % first, check whether obj is degenerate
    if isempty(V)
        E = ellipsoid;
        return;
    end
    n = size(V,1);
    [T,S,~] = svd(V);
    n_nd = rank(ellipsoid(S(:,1:n)));
    V = T'*V;
    if n_nd<n
        % remove zeros if degenerate
        V(n_nd+1:end,:) = [];
    end
    obj = mptPolytope(V');
    A = obj.P.A;
    d = obj.P.b;
    % normalize to prevent numerical issues
    fac = 1./sqrt(sum(A.^2,2));
    A = 1./sqrt(sum(A.^2,2)).*A;
    d = fac.*d;
    M = size(A,1);

    if isSolverInstalled('mosek')
       % model using mosek
        [~, res] = mosekopt('symbcon echo(0)');
        prob.bardim = 2*n_nd;
        
        gsum = @(x) x.*(x+1)./2;
    
        % number of variables for n_nd symmetric matrix
        nS = gsum(n_nd);
        % indices for lower triangular part of n_nd symmetric matrix
        [rl_S,cl_S] = find(tril(ones(n_nd))); 
        [ru_S,cu_S] = find(triu(ones(n_nd)));
        indl_S = sub2ind([n_nd,n_nd],rl_S,cl_S);
    
        % generate ids to identify same constraints: For our problem
        id_max = 0;
        id_b = id_max+(1:n_nd); id_max = max(id_b);
        id_B = id_max+(1:nS); id_max = max(id_B);
        id_Zt = id_max+(1:nS); id_max = max(id_Zt);
        % ids for diagonal entries (upper triangular matrix because 
        % transposed)
        ii_diag_u = gsum(1:n_nd);
        id_Zd = id_Zt(ii_diag_u);

        id_B_full = zeros(n_nd);
        id_B_full(indl_S) = id_B;
        id_B_full = id_B_full + tril(id_B_full,-1)';
        
        % constraints due to geomean(B) reformulation
        % formulate [B,Z;Z',diag(Z)]>=0
        % Cb
        % nothing to do...
    
        % Ab
        % entries for B
        a_subi_B_det = id_B;
        a_subj_B_det = ones(1,nS);
        a_subk_B_det = rl_S';
        a_subl_B_det = cl_S';
        a_val_B_det = -ones(1,nS);
    
        % entries for LT Z & diag(Z)
        a_subi_Z_det = [id_Zt,id_Zd];
        a_subj_Z_det = ones(1,nS+n_nd);
        a_subk_Z_det = [n_nd+ru_S',n_nd+(1:n_nd)];
        a_subl_Z_det = [cu_S',n_nd+(1:n_nd)];
        a_val_Z_det = [-ones(1,nS),-ones(1,n_nd)];
    
        % collect
        prob.bara.subi = [a_subi_B_det,a_subi_Z_det];
        prob.bara.subj = [a_subj_B_det,a_subj_Z_det];
        prob.bara.subk = [a_subk_B_det,a_subk_Z_det];
        prob.bara.subl = [a_subl_B_det,a_subl_Z_det];
        prob.bara.val = [a_val_B_det,a_val_Z_det];

        prob.barc.subj = zeros(1,0);
        prob.barc.subk = zeros(1,0);
        prob.barc.subl = zeros(1,0);
        prob.barc.val = zeros(1,0);
        
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
        id_y = [id_b,id_B,id_Zt,id_t,id_z];
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
                                    [0.5,0.5,0;0.5,-0.5,0;0,0,1];
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
                                    [0.5,0.5,0;0.5,-0.5,0;0,0,1];
            ii_t_prev = ii_t_prev + 2;
        end

        % add equalities to model ||B*A(i,:)'||_2 + A(i,:)*b <= d_i
        % or equivalently ||B*A(i,:)'||_2 <= d_i - A(i,:)*b
        n_s_p = M*(1+n_nd);
        At_p = zeros(n_s_p,size(At,2));
        % construct logical mask
        iMask = false(n_nd,size(At,2));
        iMask_b = ismember(id_y,id_b);
        for i=1:n_nd
            % id_B: mask for B*A(i,:)'
            iMask(i,ismember(id_y,id_B)) = ismember(id_B,id_B_full(i,:));
        end
        for i=1:M
            at_p_i = zeros(n_nd,size(At,2));
            for j=1:n_nd
                % set corresponding coefficients of B*A(i,:)'
                at_p_i(j,iMask(j,:)) = A(i,:);
            end
            % coefficients for d_i-A(i,:)*b
            % d_i is taken care of by prob.c (see below)
            % (set s_{..1} = d_i-A(i,:)*b & s_{...1+i} = B(i,:)*A(i,:)')
            % (i \in {1,...,n_nd})
            at_aib = zeros(1,size(At,2));
            at_aib(iMask_b) = -A(i,:);
            At_p((i-1)*(1+n_nd)+1:i*(1+n_nd),:) = [at_aib;at_p_i];
        end
        prob.a = sparse(-[At;At_p]');
        
        % objective in scalar variables:
        % we have M*(n_nd+1) rows in At (thus the same number of new scalar
        % primal variables); first of each block (of M) is supposed to be
        % d_i-A(i,:)*b, second B*A(i,:)
        % here, we only need to take care of the constant part, i.e. d_i
        % for the first block, and zeros for the remaining n_nd rows
        prob.c = [zeros(1,size(At,1)),reshape([d,zeros(M,n_nd)]',1,[])];
    
        % formulate cones
        prob.f = speye(size(prob.a,2));
        prob.g = zeros(size(prob.f,1),1);

        
        % construct 3d cones and M (n_nd+1)-dim quadratic cones
        prob.cones = [repmat([res.symbcon.MSK_CT_QUAD 3],1,n_cones),...
                    repmat([res.symbcon.MSK_CT_QUAD n_nd+1],1,M)];
    
        % inequality constraints (equality constraints just have same lb & ub)
        prob.blc = [zeros(1,n_y-1),1];
        prob.buc = prob.blc;
    
        % optimize
        [~,res] = mosekopt('minimize echo(0)',prob);

        % check if everything worked out
        if res.rcode~=0 ||  ~strcmp(res.sol.itr.solsta,'OPTIMAL')
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
        ind_b = ismember(id_y,id_b);
        b_sol = y_sol(ind_b);
        Qt = B_sol^2;
        qt = b_sol;

    elseif isSolverInstalled('sdpt3')
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
        id_b = id_max+(1:n_nd); id_max = max(id_b);
        id_B = id_max+(1:nS);
        id_B_full = zeros(n_nd);
        id_B_full(indu_S) = id_B;
        id_B_full = id_B_full + triu(id_B_full,1)';
    
        blk = cell(M+1,2);
        % blocks for M quadratic cone constraints
        blk(1:M,1) = {'q'}; blk(1:M,2) = {n_nd+1};
        % block for psd matrix (B)
        blk(M+1,:) = {'s',n_nd};
        
        C = cell(M+1,1);
        
        % y = [b;Bu_vec]
        id_y = [id_b,id_B];
        n_y = length(id_y);

        % construct quadratic cones' At

        % construct logical mask
        iMask = false(n_nd,size(n_y,2));
        for i=1:n_nd
            iMask(i,ismember(id_y,id_B)) = ismember(id_B,id_B_full(i,:));
        end

        At = cell(M+1,1);
        iMask_b = ismember(id_y,id_b);
        for i=1:M
            % constant matrix
            C{i} = sparse(1,1,d(i),n_nd+1,1);
            at_p_i = zeros(n_nd,n_y);
            for j=1:n_nd
                at_p_i(j,iMask(j,:)) = A(i,:);
            end
            % coefficients for d_i-A(i,:)*b
            % d_i is taken care of by prob.c (see below)
            % (set s_{..1} = d_i-A(i,:)*b & s_{...1+i} = B(i,:)*A(i,:)')
            % (i \in {1,...,n_nd})
            at_aib = zeros(1,n_y);
            at_aib(iMask_b) = -A(i,:);
            At{i} = [at_aib;at_p_i];
        end
        
        % objective matrix for B>=0
        C{M+1} = sparse(n_nd,n_nd);
    
        % constraint matrix for B>=0 (ensuring that "psd-B" and "y-B" are
        % equal)
        aend_subk_B = 1:nS;
        aend_subl_B = id_B;
        % extract diagonal and non-diagonal entries to scale non-diagonal
        % entries with sqrt(2)
        ind_diag_B = ismember(id_B,diag(id_B_full)');
        aend_val_B(ind_diag_B) = -ones(1,sum(ind_diag_B));
        aend_val_B(~ind_diag_B) = -sqrt(2)*ones(1,sum(~ind_diag_B));
        At{M+1} = sparse(aend_subk_B,aend_subl_B,aend_val_B,nS,n_y);
    
        % set logdet terms
        OPTIONS.parbarrier = [repmat({0},M,1);{1}];
        OPTIONS.printlevel = 0;
    
        % since sdpt3 supports logdet explicitly, b = 0
        b = sparse(n_y,1);
    
        % call solver
        [~,~,y,~,info] = sqlp(blk,At,C,b,OPTIONS);
        
        if info.termcode~=0
            throw(CORAerror('CORA:solverIssue','sdpt3'));
        end
    
        % extract solution
        B_sol_vec = y(ismember(id_y,id_B));
        B_sol = zeros(n_nd);
        B_sol(triu(ones(n_nd))==1) = B_sol_vec;
        B_sol = B_sol + triu(B_sol,1)';
        b_sol = y(ismember(id_y,id_b));
        
        Qt = B_sol^2;
        qt = b_sol;

    elseif isYalmipInstalled()
        c = sdpvar(n_nd,1);
        B = sdpvar(n_nd);
        F = cone([d-A*c,A*B]');
        options = sdpsettings;
        options.verbose = 0;
        if exist('sdpt3','file')
            options.solver = 'sdpt3';
        end
        %solve optimization problem
        diagnostics = optimize(F, -logdet(B), options);
        if diagnostics.problem ~= 0
             throw(CORAerror('CORA:solverIssue'));
        end
        Qt = value(B)^2;
        qt = value(c);
    else
        throw(CORAerror('CORA:noSuitableSolver','SDP'));
    end
    E = T*ellipsoid(blkdiag(Qt,zeros(n-n_nd)),[qt;zeros(n-n_nd,1)]);
end
% add back center
E = E + m;

%------------- END OF CODE --------------