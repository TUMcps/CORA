function E = enclosePoints(points,varargin)
% enclosePoints - enclose a point cloud with an ellipsoid
%
% Syntax:
%    E = enclosePoints(points)
%    E = enclosePoints(points,method);
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%    method - (optional) method to compute the enclosing ellipsoid
%               'cov' or 'min-vol'
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    mu = [2 3];
%    sigma = [1 1.5; 1.5 3];
%    points = mvnrnd(mu,sigma,100)';
%
%    E1 = ellipsoid.enclosePoints(points);
%    E2 = ellipsoid.enclosePoints(points,'min-vol');
%    
%    figure; hold on
%    plot(points(1,:),points(2,:),'.k');
%    plot(E1,[1,2],'r');
%    plot(E2,[1,2],'b');
%
% References:
%   [1]: Boyd; Vandenberghe: Convex Optimization
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclosePoints, interval/enclosePoints

% Authors:       Niklas Kochdumper, Victor Gassmann
% Written:       05-May-2020
% Last update:   15-March-2021
%                09-June-2022 (VG, removed false assumption on center in 'min-vol')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'cov'},varargin);

% check input arguments
inputArgsCheck({{points,'att','numeric','nonempty'};
                {method,'str',{'cov','min-vol'}}});

% remove bias
c = mean(points,2);
points = points-c;
[n,M] = size(points);

% handle degenerate case
[T,S,~] = svd(points);
n_nd = rank(S);
points = T'*points;
n_d = n-n_nd;
%r_d = zeros(n_d,1);
if n_nd<n
    % essentially, this is pca; since we know that the data is
    % "degenerate", the last row has to be ~zero, as the variance is ~0
    % remove zeros 
    % save rest to add later (don't just ignore the very small parts)
    %r_d = 1/2*(max(points(n_nd+1:end,:),[],2)-min(points(n_nd+1:end,:),[],2));
    points(n_nd+1:end,:) = [];
end

% handle special cases n_nd=0 and n_nd=1
if n_nd==0
    % all zero matrix
    E = ellipsoid.empty(n);

elseif n_nd==1
    % interval arithmetic (is exact in this case)
    r = 1/2*(max(points)-min(points));
    q = 1/2*(max(points)+min(points));
    E = ellipsoid(r^2,q);
    
elseif strcmp(method,'cov')

    % compute the covariance matrix
    C = cov(points');

    % singular value decomposition
    [U,~,~] = svd(C);
    
    % required ellipsoid radius in transformed space
    orientedMatrix=U'*points;
    m1 = max(orientedMatrix,[],2);
    m2 = min(orientedMatrix,[],2);

    nt = max([m1,m2],[],2);
    
    % enclosing ellipsoid in the transformed space
    B = diag(nt.^2);
    
    maxDist = max(sum(diag(1./nt.^2)*orientedMatrix.^2,1));
    B = B * maxDist;
    
    E = ellipsoid(B);
    
    % transform back to original space
    E = U*E;
    
elseif strcmp(method,'min-vol')

    persistent isMosek
    if isempty(isMosek)
        isMosek = isSolverInstalled('mosek');
    end
    persistent isSDPT3
    if isempty(isSDPT3)
        isSDPT3 = isSolverInstalled('sdpt3');
    end

    if isMosek
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

        % add equalities to model B*x^{(i)}+b = s_i
        % we need an additional "1" to later model the quadratic cone for
        % each point, thus for each point in "points", we need 1+n rows
        % e.g. for B = [B11,B12;B12,B22], and x^{(i)} = [x1;x2], we need the
        % constraint: 
        %%% [1;0;0] + [0;b1;b2] + [0;B11*x1+B12*x2;B12*x1+B22*x2] = s_p_i,
        % where s_p_i are cone variables introduced later. Thus
        n_s_p = M*(1+n_nd);
        At_p = zeros(n_s_p,size(At,2));
        % construct logical mask
        iMask = false(n_nd,size(At,2));
        for i=1:n_nd
            % id_B
            % i-th row has B(i,:)*x^{(i)} and b_i
            iMask(i,ismember(id_y,id_b(i))) = true;
            iMask(i,ismember(id_y,id_B)) = ismember(id_B,id_B_full(i,:));
        end
        for i=1:M
            at_p_i = zeros(n_nd,size(At,2));
            for j=1:n_nd
                at_p_i(j,iMask(j,:)) = [-1,-points(:,i)'];
            end
            % here we set s_{..1} = 1 & s_{...1+i} = B(i,:)*x^{(i)}+b_i
            % (i \in {1,...,n_nd})
            At_p((i-1)*(1+n_nd)+1:i*(1+n_nd),:) = [zeros(1,size(At,2));at_p_i];
        end
        prob.a = sparse(-[At;At_p]');
        
        % objective in scalar variables
        prob.c = [zeros(1,size(At,1)),repmat([1,zeros(1,n_nd)],1,M)];
    
        % formulate cones
        prob.f = speye(size(prob.a,2));
        prob.g = zeros(size(prob.f,1),1);

        
        % construct 3d cones and M (n_nd+1)-dim quadratic cones
        prob.cones = [repmat([res.symbcon.MSK_CT_QUAD 3],1,n_cones),...
                    repmat([res.symbcon.MSK_CT_QUAD n_nd+1],1,M)];
    
        % inequality constraints (equality constraints just have same lb and ub)
        prob.blc = [zeros(1,n_y-1),1];
        prob.buc = prob.blc;
    
        % optimize
        [~,res] = mosekopt('minimize echo(0)',prob);

        % check if everything worked out
        if res.rcode~=0 || ~strcmp(res.sol.itr.solsta,'OPTIMAL')
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
        E = ellipsoid(inv(B_sol^2),-B_sol\b_sol);

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
            iMask(i,ismember(id_y,id_b(i))) = true;
            iMask(i,ismember(id_y,id_B)) = ismember(id_B,id_B_full(i,:));
        end

        At = cell(M+1,1);
        for i=1:M
            % constant matrix
            C{i} = sparse(1,1,1,n_nd+1,1);
            at_p_i = zeros(n_nd,n_y);
            for j=1:n_nd
                at_p_i(j,iMask(j,:)) = [-1,-points(:,i)'];
            end
            At{i} = [zeros(1,n_y);at_p_i];
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
        E = ellipsoid(inv(B_sol^2),-B_sol\b_sol);

    elseif isYalmipInstalled()
        nt = size(points,1);
        B = sdpvar(nt);
        b = sdpvar(nt,1);
        X = points;
        F = [];
        for i=1:M
            F = [F,norm(B*points(:,i)+b)<=1];
        end
        options = sdpsettings;
        options.verbose = 0;
        % solve optimization problem
        diagnostics = optimize(F, -geomean(B), options);
        warning("YALMIP was used to model the problem - " + ...
            "consider installing a supported solver to speed up computation...");
        if diagnostics.problem ~= 0
            throw(CORAerror('CORA:solverIssue'));
        end
        E = ellipsoid(inv(value(B)^2),-value(B)\value(b));
        
    else
        throw(CORAerror('CORA:noSuitableSolver','SDP'));
    end
end

% backtransform and add back center
E_ext = E;
if n_d>0
    E_ext = ellipsoid(blkdiag(E_ext.Q,zeros(n_d)),[E_ext.q;zeros(n_d,1)]);
end
E = T*E_ext + c;

% ------------------------------ END OF CODE ------------------------------
