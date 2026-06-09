function Z = reduceUnderApprox(Z,method,order)
% reduceUnderApprox - reduces the order of a zonotope so that an
%    under-approximation of the original set is obtained
%
% Syntax:
%    Z = reduceUnderApprox(Z,method,order)
%
% Inputs:
%    Z - zonotope object
%    method - reduction method ('sum','scale','linProg','wetzlinger')
%    order - zonotope order
%
% Outputs:
%    Z - reduced zonotope
%
% Example: 
%    Z = zonotope([1;-1],[3 2 -3 -1 2 4 -3 -2 1; 2 0 -2 -1 2 -2 1 0 -1]);
%
%    Zsum = reduceUnderApprox(Z,'sum',3); 
%    Zscale = reduceUnderApprox(Z,'scale',3);
%    ZboxLP = reduceUnderApprox(Z,'linProg',3);
%    ZboxWetz = reduceUnderApprox(Z,'wetzlinger',3);
%   
%    figure; hold on;
%    plot(Z,[1,2],'r','LineWidth',2);
%    plot(Zsum,[1,2],'b');
%    plot(Zscale,[1,2],'g');
%    plot(ZboxLP,[1,2],'m');
%    plot(ZboxWetz,[1,2],'c');
%
% References:
%    [1] Sadraddini et al. "Linear Encodings for Polytope Containment
%        Problems", CDC 2019
%    [2] Wetzlinger et al. "Adaptive Parameter Tuning for Reachability 
%        Analysis of Nonlinear Systems", HSCC 2021  
%    [3] Alizadeh and Goldfarb. "Second-Order Cone Programming",
%        Mathematical Programming 2003
%    [4] Yang and Ozay, "Scalable zonotopic under-approximation of 
%        backward reachable sets for uncertain linear systems," IEEE 
%        Control Syst. Lett., vol. 6, 2022. 
%    [5] Raghuraman and Koeln, "Set operations and order reductions for 
%        constrained zonotopes," Automatica, vol. 139, 2022. 
%    [6] Kochdumper and Bak, "Conformant synthesis for Koopman operator 
%        linearized control systems," in CDC 2022
%    [7] Luetzow et al., "Underapproximative Methods for the Order 
%        Reduction of Zonotopes," IEEE Control Syst. Lett., vol. 9, 2025. 
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: reduce

% Authors:       Niklas Kochdumper
% Written:       19-November-2018
% Last update:   29-August-2019
%                15-April-2020 (added additional reduction techniques)
%                16-January-2025 (LL, NK, rename and add new methods)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    inputArgsCheck({{Z,'att','zonotope'};
                    {method,'str',{'sadraddini','yang','raghuraman','kochdumper',...
                    'scale','boxLP','boxWetzlinger','boxCone','nlp','cluster', ...
                    'sum','linProg','wetzlinger'}}; % old names
                    {order,'att','numeric','nonnan'}});
    
    % remove all-zero generators
    Z = compact_(Z,'zeros',eps);

    % check if reduction is required
    [n, nrOfGens] = size(generators(Z));

    if n*order < nrOfGens
        
        % reduce with the selected method
        % methods from literature
        if strcmp(method,'sadraddini') % according to [1]
            Z = aux_reduceUnderApproxSadraddini(Z,order);
        elseif strcmp(method,'yang') % according to [4]
            Z = aux_reduceUnderApproxYang(Z,order);
        elseif strcmp(method,'raghuraman') % according to [5]
            Z = aux_reduceUnderApproxRaghuraman(Z,order);
        elseif strcmp(method,'kochdumper') ...  % according to [6]
                || strcmp(method,'sum') % old name
            Z = aux_reduceUnderApproxKochdumper(Z,order);
        % scale method
        elseif strcmp(method,'scale') % according to [7]
            Z = aux_reduceUnderApproxScale(Z,order);
        % box methods
        elseif strcmp(method,'boxCone') % according to [7]
            Z = aux_reduceUnderApproxBoxCone(Z,order);
        % NLP methods
        elseif strcmp(method,'nlp') % according to [7]
            Z = aux_reduceUnderApproxNLP(Z,order);
        % clustering methods
        elseif strcmp(method,'cluster') % according to [7]
            Z = aux_reduceUnderApproxCluster(Z,order);
        end
        % LEGACY methods 
        elseif strcmp(method,'boxLP') ... % new name
                || strcmp(method,'linProg') % old name
            Z = aux_reduceUnderApproxBoxLP(Z,order);
        elseif strcmp(method,'boxWetzlinger') ... % new name
                || strcmp(method,'wetzlinger') % old name
            Z = aux_reduceUnderApproxBoxWetzlinger(Z,order);
        elseif strcmp(method,'boxNLP')
            Z = aux_reduceUnderApproxBoxNLP(Z,order);
    else
        return;
    end

end


% Auxiliary functions -----------------------------------------------------

function Zred = aux_reduceUnderApproxSadraddini(Z,order)
% reduction method in Corrollary 6 in [1]

    c = Z.c; G = Z.G;
    [n,m] = size(G); mred = floor(n*order);

    % randomly initialize the generator matirx of the reduced zonotope
    Gred = 2*(rand(n,mred)-0.5);

    % compute feasible initial generator matrix by scaling the generator
    % matirx accordingly
    [A,b,Aeq,beq,lb,ub,f] = aux_contConstrTedrake(Gred,G);

    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    problem.options = options;
    
    [x,~,exitflag,output] = CORAlinprog(problem);
    
    if exitflag <= 0
        fprintf("Sadraddini: " + output.message + "\n");
        Zred = zonotope(c,[]);
        return 
        %throw(CORAerror('CORA:solverIssue'));
    end
    
    lambda = x(1:mred);
    Gred = Gred*diag(lambda);

    % iteratively improve the solution by alternatingly solving for the
    % optimal Gred and T1
    for i = 1:10

        % solve Equation (43) in [1] for the optimal T1
        [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniFixGred(G,Gred);

        problem.f = f;
        problem.Aineq = A; problem.bineq = b;
        problem.Aeq = Aeq; problem.beq = beq;
        problem.lb = lb; problem.ub = ub;
        
        [x,~,exitflag,output] = CORAlinprog(problem);
        
        if exitflag <= 0
            fprintf("Sadraddini: " + output.message + "\n");
            Gred = [];
            break
            %throw(CORAerror('CORA:solverIssue'));
        end

        ind = 1 + m*mred + (1:m*mred);
        T1 = reshape(x(ind),[mred,m]);

        % solve Equation (43) in [1] for the optimal Gred
        [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniFixT1(G,T1);

        problem.f = f;
        problem.Aineq = A; problem.bineq = b;
        problem.Aeq = Aeq; problem.beq = beq;
        problem.lb = lb; problem.ub = ub;
        
        [x,~,exitflag, output] = CORAlinprog(problem);
        
        if exitflag <= 0
            fprintf("Sadraddini: " + output.message + "\n");
            Gred = [];
            break
            %throw(CORAerror('CORA:solverIssue'));
        end

        ind = 1 + (1:n*mred);
        Gred = reshape(x(ind),[n,mred]);
    end

    % construct reduced order zonotope
    Zred = zonotope(c,Gred);
end

function Zred = aux_reduceUnderApproxYang(Z,order)
% iteratively find and combine closely aligned generators

G = generators(Z);
closeness_ij = [];
while size(G,2) > order*dim(Z)
    eta = size(G,2);

    % find closely aligned generators g_i and g_j
    if isempty(closeness_ij)
        % compute matrix evaluating the closeness of g_i and g_j from scratch
        closeness_ij = inf*ones(eta,eta);
        for i=1:eta-1
            g_i = G(:,i);
            for j=i+1:eta
                g_j = G(:,j);
                closeness_ij(i,j) = norm(g_i) * norm(g_j - 1/(norm(g_i)^2)*g_i*g_j'*g_i );
            end
        end
        %eval_norm2 = sqrt(sum(G'.^2,2)). * sqrt(sum(G'- 1/sum(G'.^2,2).* ???,2))
    else
        % remove i'th and j'th rows and column in norm_ij
        closeness_ij(:,[i, j]) = [];
        closeness_ij([i, j],:) = [];

        % only evaluate closeness for the newly added generator
        closeness_ij = [closeness_ij inf*ones(size(closeness_ij,1),1); ...
            inf*ones(1,size(closeness_ij,2)+1)]; % add new row and column
        j = eta;
        g_j = G(:,j); 
        for i=1:eta-1
            g_i = G(:,i);
            closeness_ij(i,j) = norm(g_i) * norm(g_j - 1/(norm(g_i)^2)*g_i*g_j'*g_i );
        end
    end

    [~,idx] = min(closeness_ij, [], 'all');
    [i,j] = ind2sub(size(closeness_ij),idx);
    g_i = G(:,i);
    g_j = G(:,j);

    % remove g_i and g_j and add combined vector 
    G(:,[i, j]) = [];
    G_r = G'/(G*G');
    if norm((g_i+g_j)'* G_r')/norm((g_i-g_j)'* G_r') >= 1
        g_new = g_i+g_j;
    else
        g_new = g_i-g_j;
    end
    G = [G g_new];
end

Zred = zonotope([center(Z),G]);
end

function Zred = aux_reduceUnderApproxRaghuraman(Z,order)
% add or subtract small generators to/from mostly aligned big generators

    % select generators to reduce
    n = dim(Z);
    eta_new = floor(order*n);
    
    [c,G_big,G_small] = aux_selectSmallestGenerators(Z,eta_new,'length');

    eta = size(generators(Z),2);
    G_prod = G_small'*G_big;
    alphas = abs(G_prod);
    [alpha_max,idx] = max(alphas,[],2);
    T = zeros(eta-eta_new, eta_new);
    
    T(sub2ind(size(T), 1:eta-eta_new, idx')) = ...
        G_prod(sub2ind(size(G_prod), 1:eta-eta_new, idx'))./alpha_max';
    T = [eye(eta_new); T];

    % construct the reduced zonotope object
    Zred = zonotope([c,generators(Z)*T]);
end

function Zred = aux_reduceUnderApproxKochdumper(Z,order)
% sum up the generators that are reduced to obtain an inner-approximation

    % select generators to reduce
    n = dim(Z);
    N = floor(order*n - 1);
   
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'length');

    % under-approximate the generators that are reduced by one generator
    % corresponding to the sum of generators
    g = sum(Gred*diag(sign(Gred(1,:))),2);

    % construct the reduced zonotope object
    Zred = zonotope([c,G,g]);
end

function Zred = aux_reduceUnderApproxScale(Z,order)
% An over-approximative reduced zonotope is computed first. This zonotope
% is then scaled using linear programming until it is fully contained in 
% the original zonotope 

    % over-approximative reduction of the zonotope
    Zo = reduce(Z,'constOpt',order);
    
    % get conditions for linear program to scale the over-approximative 
    % zonotope until it is contained inside the original zonotope
    [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrTedrake(Zo.G,Z.G);

    % init linprog struct
    problem.f = f';
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    % solve linear program
    [x,~,exitflag,output] = CORAlinprog(problem);
    
    if exitflag <= 0
        fprintf("Scale: " + output.message + "\n");
        Gred = [];
        %throw(CORAerror('CORA:solverIssue'));
    else
        lambda = x(ind);
        Gred = Zo.G*diag(lambda);
    end

    
    % construct final zonotope
    Zred = zonotope([Z.c,Gred]);
end

function Zred = aux_reduceUnderApproxBoxLP(Z,order)
% reduce the zonotope order by computing an interval under-approximation of
% the zonotope spanned by the reduced generators using linear programming

    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'box');

    % principal component analysis (PCA)
    p = [-Gred,Gred];
    C = cov(p');
    [U,~,~] = svd(C);
    Gred = U'*Gred;

    % obtain constraints for zonotope in zonotope containment according to
    % Theorem 3 in [2]
    [A,b,Aeq,beq,lb,ub,f] = aux_contConstrTedrake(eye(n),Gred);

    % init linprog struct
    problem.f = f;
    problem.Aineq = A; problem.bineq = b;
    problem.Aeq = Aeq; problem.beq = beq;
    problem.lb = lb; problem.ub = ub;
    
    % solve linear program
    [x,~,exitflag,output] = CORAlinprog(problem);

    if exitflag <= 0
        fprintf("BoxLP: " + output.message + "\n");
        Zred = zonotope([c,G]);
        %throw(CORAerror('CORA:solverIssue'));
    else
        gamma = x(1:n);

        % construct the reduced zonotope
        Zred = zonotope([c,G,U*diag(gamma)]);
    end
end

function Zred = aux_reduceUnderApproxBoxWetzlinger(Z,order)
% reduction based on the Hausdorff distance between a zonotope and its
% interval enclosure (see Theorem 3.2 in [2])

    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'box');

    % principal component analysis (PCA)
    p = [-Gred,Gred];
    C = cov(p');
    [U,~,~] = svd(C);
    Gred = U'*Gred;

    % compute over-approximation of the Hausdorff distance between the
    % zonotope and its box enclsoure according to Theorem 3.2 in [2]
    Gabs = abs(Gred);
    G_hat = Gabs;
    
    for k = 1:size(G_hat,2)
        [~,i_star] = max(Gabs(:,k));
        G_hat(i_star(1),k) = 0;
    end
    
    e = 2*sum(G_hat,2);
    
    % subtract error bound from the box enclosure
    gamma_o = sum(Gabs,2);
    
    if any(gamma_o - e <= 0)
        gamma = zeros(n,1);
    else
        gamma = gamma_o - e;
    end
    
    % combine the box inner-approximation with the unreduced generators 
    Zred = zonotope(c,[G,U*diag(gamma)]);
end

function Zred = aux_reduceUnderApproxBoxNLP(Z,order)
% reduce the zonotope order by computing an interval under-approximation of
% the zonotope spanned by the reduced generators using nonlinear programming
   
    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'box');

    % principal component analysis (PCA)
    p = [-Gred,Gred];
    C = cov(p');
    [U,~,~] = svd(C);
    Gred = U'*Gred;

    % obtain constraints for zonotope in zonotope containment according to
    % Theorem 3 in [2]
    [A,b,Aeq,beq,lb,ub] = aux_contConstrTedrake(eye(n),Gred);

    % cost function and initial guess
    cost_fun = @(x) aux_cost_volumeInterval(x,n);
    x0 = [0.1*sum(abs(Gred),2);zeros(size(A,2)-n,1)];

    % solve nonlinear program
    persistent options
    if isempty(options)
        options = optimoptions('fmincon','Algorithm','sqp', ...
            'MaxIterations',100000,'display','off',...
            'MaxFunctionEvaluations',Inf,'SpecifyObjectiveGradient',true);
    end
    [x, ~, exitflag,output] = fmincon(cost_fun,x0,A,b,Aeq,beq,lb,ub,[],options);
    
    if exitflag <= 0
        fprintf("BoxNLP: " + output.message + "\n");
        Zred = zonotope([c,G]);
        %throw(CORAerror('CORA:solverIssue'));
    else
        gamma = x(1:n);

        % construct the reduced zonotope
        Zred = zonotope([c,G,U*diag(gamma)]);
    end
end

function Zred = aux_reduceUnderApproxBoxCone(Z,order)
% reduce the zonotope order by computing an interval under-approximation of
% the zonotope spanned by the reduced generators using second-order cone 
% programming
   
    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'box');

    % principal component analysis (PCA)
    p = [-Gred,Gred];
    C = cov(p');
    [U,~,~] = svd(C);
    Gred = U'*Gred;

    % obtain constraints for zonotope in zonotope containment according to
    % Theorem 3 in [2]
    [A,b,Aeq,beq,lb,ub,~] = aux_contConstrTedrake(eye(n),Gred);

    % second-order cone constraints for optmizing the volume according to
    % Section 2.3.e) in [3]
    ind = 1:n; len = size(A,2);

    [socCon,len_,Aeq_,beq_] = aux_constraintsGeoMean(ind,len);

    f = zeros(1,len_); f(end) = -1;
    A = blkdiag(A,-eye(len_-len)); b = [b;zeros(len_-len,1)];
    Aeq = [[Aeq,zeros(size(Aeq,1),len_-len)]; Aeq_]; beq = [beq;beq_];
    
    % solve second-order cone program
    persistent options
    if isempty(options)
        options = optimoptions('coneprog','Display','off');
    end

    [x,~,exitflag,output] = coneprog(f,socCon,A,b,Aeq,beq,lb,ub,options);

    if exitflag <= 0
        fprintf("BoxCone: " + output.message + "\n");
        Zred = zonotope([c,G]);
        %throw(CORAerror('CORA:solverIssue'));
    else
        gamma = x(1:n);

        % construct the reduced zonotope
        Zred = zonotope([c,G,U*diag(gamma)]);
    end
end

function [c,g] = aux_cost_volumeInterval(x,n)
% cost function for maximizing the volume of an interval together with the
% corresponding gradient

    % cost function
    s = x(1:n);
    c = -prod(s);

    % gradient
    g = zeros(size(x));

    for i = 1:n
        ind = setdiff(1:n,i);
        g(i) = -prod(s(ind));
    end

    if all(g(1:n) == 0)
        g(1:n) = -ones(n,1);
    end
end

function Zred = aux_reduceUnderApproxNLP(Z,order)
% reduce the zonotope order by computing an interval under-approximation of
% the zonotope spanned by the reduced generators using linear programming

    % select generators to reduce
    n = dim(Z);

    % number of generators which are copied directly
    orderOpt = 1;
    N = floor(order*n-orderOpt*n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N,'length');
    eta = size(Gred,2);
    eta_r = order*n - size(G,2);

    % constraint || T ||_inf <= 1
    M = reshape(1:eta*eta_r,[eta, eta_r]); 
    [A,b,Aeq,beq] = aux_infinityNormConstraint(M,eta*eta_r);

    % constraints -1<=T_ij<=1, 0<=mu<=1 
    % (A,b encode 0<=mu -> not required anymore)
    lb = [-1* ones(eta*eta_r,1); zeros(2*eta*eta_r,1)];
    ub = 1* ones(3*eta*eta_r,1);

    % initial guess: copy the first eta_r generators of G_red, add all
    % other generators to the last copied generator
    T0 = [eye(eta_r); [zeros(eta-eta_r, eta_r-1) ones(eta-eta_r,1)]];
    Mu0 = [T0 zeros(size(T0))];
    x0 = [reshape(T0,[eta*eta_r,1]); reshape(Mu0',[2*eta*eta_r,1]);];

    GredSparseCells = repmat({sparse(Gred)},[eta_r,1]);
    GredBig = blkdiag(GredSparseCells{:});
    cost_fun = @(x) aux_cost_volume(x,Gred,eta,eta_r,GredBig);
    
    % solve nonlinear program
    persistent options
    if isempty(options)
        options = optimoptions('fmincon','Algorithm','sqp', ...
            'MaxIterations',100000,'display','off',...
            'MaxFunctionEvaluations',Inf,'SpecifyObjectiveGradient',true);
    end
    [x, ~, exitflag,output] = fmincon(cost_fun,x0,[],[],Aeq,beq,lb,ub,[],options);
    if exitflag <= 0
        fprintf("NLP: " + output.message + "\n");        
        x = x0;
    end

    if cost_fun(x) > cost_fun(x0)
        x = x0;
        fprintf("Reset optimization result to x0 \n");
    end
    
    T = reshape(x(1:eta*eta_r), [eta, eta_r]);

    % construct the reduced zonotope
    Zred = zonotope([c,G,Gred*T]);
end

function [c,g] = aux_cost_volume(x,Gred,eta,eta_r,GredBig)
% cost function for volume optimization together with corresponding gradient

    % compute cost
    Gamma = reshape(x(1:eta*eta_r), [eta, eta_r]);
    d = det(Gred*Gamma);

    c = -abs(d);

    % compute gradient
    g = zeros(size(x));
    adjointGrad = sum(reshape(adjoint(Gred*Gamma)',[eta_r^2,1]).*GredBig,1);
    g(1:eta*eta_r) = -sign(d)*adjointGrad;
end

function Zred = aux_reduceUnderApproxCluster(Z,order)
% cluster generators and sum similar ones

G = generators(Z);

norm_factor = G(1,:);
norm_factor(norm_factor==0)=1;
G_norm = G ./ norm_factor; % normalize by the first element of each generator
idx = kmeans(G_norm', order*dim(Z), 'Distance', 'cosine');
G_new = zeros(dim(Z), order*dim(Z));
for i = 1:length(idx)
    i_cluster = idx(i);
    G_new(:,i_cluster) = G_new(:,i_cluster) + G(1,i)/abs(norm_factor(1,i)) * G(:,i);
end

Zred = zonotope(center(Z), G_new);
end

function [c,G,Gred] = aux_selectSmallestGenerators(Z,N,type)
% select the generators that are reduced

    % obtain object properties
    c = Z.c;
    G_ = Z.G;

    % sort according to generator length
    if strcmp(type,'length')
        genPriority = sum(G_.^2,1);
    elseif strcmp(type,'box')
        Gabs = abs(G_);
        genPriority = sum(Gabs,1) - max(Gabs,[],1);
    end

    [~,ind] = sort(genPriority,'descend');
    
    % split into reduced and unreduced generators
    G = G_(:,ind(1:N));
    Gred = G_(:,ind(N+1:end));
end


function [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrTedrake(Gx,Gy)
% this function returns sufficient conditions for the zonotope X to be
% contained in zonotope Y according to Theorem 3 in [1]. The optimization 
% variable x is defined as follows:     
%
%   x = [s,T(:,1),...,T(:,nx),mu],
%
% where s >= 0 are the scaling factors and T and mu are
% auxiliary variables

    nx = size(Gx,2); ny = size(Gy,2);

    % constraint Gx*diag(s) = Gy*T
    GyBlkCells = repmat({Gy},[nx,1]);
    A1_ = blkdiag(GyBlkCells{:});

    GxColCells = num2cell(Gx,1);
    A2_ = blkdiag(GxColCells{:});
    
    Aeq1 = [A2_,-A1_];
    beq1 = zeros(size(Aeq1,1),1);

    % constraint s > 0
    A1 = [-eye(nx),zeros(nx,size(Aeq1,2)-nx)];
    b1 = zeros(nx,1);

    % constraint || T ||_inf <= 1
    M = reshape(nx+1:nx+ny*nx,[ny,nx]); 

    [A2,b2,Aeq2,beq2] = aux_infinityNormConstraint(M,nx+ny*nx);
    
    % objective function
    f = zeros(size(A2,2),1);
    f(1:nx,1) = -sqrt(sum(Gx.^2,1))';
    
    % assemble overall constraint matrices
    beq = [beq1;beq2];
    Aeq = [[Aeq1,zeros(size(Aeq1,1),size(Aeq2,2)-size(Aeq1,2))];Aeq2];
    
    A = [[A1,zeros(size(A1,1),size(A2,2)-size(A1,2))];A2];
    b = [b1;b2];
    
    lb = []; ub = [];
    
    ind = 1:nx;
end

function [A,b,Aeq,beq] = aux_infinityNormConstraint(M,len)
% encoding of the infinity norm contraint ||M||_inf <= 1, where the
% matrix M stores the indices of the corresponding variables and len is the
% length of the vector of optimization variables

    [n,m] = size(M);
    A = []; b = []; Aeq = []; beq = [];

    % construct vertices of the L1-norm cube
    V = [eye(m),-eye(m)];

    % loop over all rows of the matrix M
    for i = 1:n

        ind = len + (((i-1)*2*m + 1):i*2*m);

        % constraint M(i,:)^T = sum_{j=1}^2m mu_j V(:,j)
        Aeq1 = zeros(m,len + 2*n*m); beq1 = zeros(m,1);

        Aeq1(:,M(i,:)) = -eye(m);
        Aeq1(:,ind) = V;

        % constraint sum_{j=1}^2m mu_j = 1
        Aeq2 = zeros(1,len + 2*n*m); beq2 = 1;

        Aeq2(:,ind) = ones(1,length(ind));

        % constraint mu_j >= 0
        A1 = zeros(2*m,len + 2*n*m); b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        % combine with previous matrices
        Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        A = [A;A1]; b = [b;b1];
    end
end

function [A,b,Aeq,beq] = aux_infinityNormConstraintVar(M,len,d)
% encoding of the infinity norm contraint ||M||_inf <= delta, where the
% matrix M stores the indices of the corresponding variables, d is the 
% index of the variable delta, and len is the length of the vector of 
% optimization variables

    [n,m] = size(M);
    A = []; b = []; Aeq = []; beq = [];

    % construct vertices of the L1-norm cube
    V = [eye(m),-eye(m)];

    % loop over all rows of the matrix M
    for i = 1:n

        ind = len + (((i-1)*2*m + 1):i*2*m);

        % constraint M(i,:)^T = sum_{j=1}^2m mu_j V(:,j)
        Aeq1 = zeros(m,len + 2*n*m); beq1 = zeros(m,1);

        Aeq1(:,M(i,:)) = -eye(m);
        Aeq1(:,ind) = V;

        % constraint sum_{j=1}^2m mu_j = delta
        Aeq2 = zeros(1,len + 2*n*m); beq2 = 0;

        Aeq2(:,d) = -1;
        Aeq2(:,ind) = ones(1,length(ind));

        % constraint mu_j >= 0
        A1 = zeros(2*m,len + 2*n*m); b1 = zeros(2*m,1);

        A1(:,ind) = -eye(2*m);

        % combine with previous matrices
        Aeq = [Aeq;Aeq1;Aeq2]; beq = [beq;beq1;beq2];
        A = [A;A1]; b = [b;b1];
    end
end

function [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniFixT1(G,T1)
% This function returns the optimization problem in Corollary 6 in [1] for 
% computing the optimal zonotope inner-approximation. The optimization 
% variables are    
%
%   x = [delta; Gred(:,1); ... ; Gred(:,mred); ...
%                  T0(:,1); ... ; T0(:,mred); Delta(:,1); ... ; Delta(:,m)]
%
% since T1 is fixed

    [n,m] = size(G);
    mred = size(T1,1);

    % constraint Gred = G*T0
    GblkCells = repmat({G},[mred,1]);

    Aeq1 = [zeros(n*mred,1),-eye(n*mred),blkdiag(GblkCells{:}),zeros(n*mred,n*m)];
    beq1 = zeros(n*mred,1);

    % constraint G = Gred*T1 + Delta
    Gred_ = reshape(1:n*mred,[n,mred]);
    AbilinT1 = zeros(n*m,n*mred); cnt = 1;

    for i = 1:m
        for j = 1:n
            AbilinT1(cnt,Gred_(j,:)) = T1(:,i);
            cnt = cnt + 1;
        end
    end

    Aeq2 = [zeros(n*m,1),AbilinT1,zeros(n*m,m*mred),eye(n*m)];
    beq2 = reshape(G,[n*m,1]);

    % constraint || T0 ||_inf \leq 1
    M = 1 + n*mred + reshape(1:m*mred,[m,mred]);

    [A1,b1,Aeq3,beq3] = aux_infinityNormConstraint(M,size(Aeq2,2));

    % constraint || \Delta ||_inf \leq \delta
    M = 1 + n*mred + m*mred + reshape(1:n*m,[n,m]);

    [A2,b2,Aeq4,beq4] = aux_infinityNormConstraintVar(M,size(Aeq3,2),1);

    % constraint delta >= 0
    A3 = zeros(1,size(A2,2)); A3(1) = -1; b3 = 0;
    
    % objective minimize delta
    f = zeros(1,size(A2,2)); f(1) = 1;

    % combine constraint matrices
    Aeq = [Aeq1;Aeq2];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq3,2)-size(Aeq,2))];Aeq3];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq4,2)-size(Aeq,2))];Aeq4];

    A = [[A1,zeros(size(A1,1),size(A2,2)-size(A1,2))];A2];
    A = [A;A3];

    beq = [beq1;beq2;beq3;beq4];
    b = [b1;b2;b3];

    lb = []; ub = [];
end

function [A,b,Aeq,beq,lb,ub,f] = aux_optProbSadraddiniFixGred(G,Gred)
% This function returns the optimization problem in Corollary 6 in [1] for 
% computing the optimal zonotope inner-approximation. The optimization 
% variables are    
%
%   x = [delta; T0(:,1); ... ; T0(:,mred); 
%               T1(:,1); ... ; T1(:,m); Delta(:,1); ... ; Delta(:,m)]
%
% since Gred is fixed.

    [n,m] = size(G);
    mred = size(Gred,2);

    % constraint Gred = G*T0
    GblkCells = repmat({G},[mred,1]);

    Aeq1 = [zeros(n*mred,1),blkdiag(GblkCells{:}),zeros(n*mred,mred*m+n*m)];
    beq1 = -reshape(Gred,[n*mred,1]);

    % constraint G = Gred*T1 + Delta
    GredBlkCells = repmat({Gred},[m,1]);

    Aeq2 = [zeros(n*m,1),zeros(n*m,mred*m),blkdiag(GredBlkCells{:}),eye(n*m)];
    beq2 = -reshape(G,[n*m,1]);

    % constraint || T0 ||_inf \leq 1
    M = 1 + reshape(1:m*mred,[m,mred]);

    [A1,b1,Aeq3,beq3] = aux_infinityNormConstraint(M,size(Aeq2,2));

    % constraint || T1 ||_inf \leq 1
    M = 1 + m*mred + reshape(1:m*mred,[mred,m]);

    [A2,b2,Aeq4,beq4] = aux_infinityNormConstraint(M,size(Aeq3,2));

    % constraint || \Delta ||_inf \leq \delta
    M = 1 + m*mred + m*mred + reshape(1:n*m,[n,m]);

    [A3,b3,Aeq5,beq5] = aux_infinityNormConstraintVar(M,size(Aeq4,2),1);

    % constaint delta >= 0
    A4 = zeros(1,size(A3,2)); A4(1) = -1; b4 = 0;
    
    % objective minimize delta
    f = zeros(1,size(A4,2)); f(1) = 1;

    % combine constraint matrices
    Aeq = [Aeq1;Aeq2];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq3,2)-size(Aeq,2))];Aeq3];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq4,2)-size(Aeq,2))];Aeq4];
    Aeq = [[Aeq,zeros(size(Aeq,1),size(Aeq5,2)-size(Aeq,2))];Aeq5];

    A = [[A1,zeros(size(A1,1),size(A2,2)-size(A1,2))];A2];
    A = [[A,zeros(size(A,1),size(A3,2)-size(A,2))];A3];
    A = [A;A4];

    beq = [beq1;beq2;beq3;beq4;beq5];
    b = [b1;b2;b3;b4];

    lb = []; ub = [];
end

function [socCon,cnt,Aeq,beq] = aux_constraintsGeoMean(ind,len)
% compute the second order cone constraints required for maximizing the
% geometric mean accoding to Section 2.3.e) in [3]

    % create additional fake scaling factors s_i to match 2^n
    aux = len + (1:(2^ceil(log2(length(ind))) - length(ind)));

    if ~isempty(aux)
        ind = [ind,aux]; len = aux(end);
    end

    % generate all constraints
    cnt = len; con = {};

    while length(ind) > 1

        indNew = [];

        % replace s_i*s_j by w^2 <= x*y
        for i = 1:floor(length(ind)/2)
            con{end+1}.w = cnt + 1;
            con{end}.x = ind(2*(i-1)+1);
            con{end}.y = ind(2*(i-1)+2);
            cnt = cnt + 1;
            indNew = [indNew,con{end}.w];
        end

        ind = indNew;
    end

    % convert costraints w^2 <= x*y to quadratic cone constraint 
    % || [2w; x-y] ||_2 <= x + y according to start of Section 2.3 in [3]
    for i = 1:length(con)

        A = zeros(2,cnt); b = zeros(2,1); d = zeros(1,cnt); g = 0;

        w = con{i}.w; x = con{i}.x; y = con{i}.y;
        A(1,w) = 2; A(2,x) = 1; A(2,y) = -1; d(x) = 1; d(y) = 1;

        socCon(i) = secondordercone(A,b,d,g);
    end

    % constrain the fake scaling factors to s_i = 1
    Aeq = zeros(length(aux),cnt); beq = ones(length(aux),1);
    Aeq(1:length(aux),aux) = eye(length(aux));
end

% ------------------------------ END OF CODE ------------------------------
