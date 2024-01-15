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
%    ZlinProg = reduceUnderApprox(Z,'linProg',3);
%   
%    figure; hold on;
%    plot(Z,[1,2],'r','LineWidth',2);
%    plot(Zsum,[1,2],'b');
%    plot(Zscale,[1,2],'g');
%    plot(ZlinProg,[1,2],'m');
%
% References:
%    [1] Sadraddini et al. "Linear Encodings for Polytope Containment
%        Problems", CDC 2019
%    [2] Wetzlinger et al. "Adaptive Parameter Tuning for Reachability 
%        Analysis of Nonlinear Systems", HSCC 2021             
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    inputArgsCheck({{Z,'att','zonotope'};
                    {method,'str',{'sum','scale','linProg','wetzlinger'}};
                    {order,'att','numeric','nonnan'}});
    
    % remove all-zero generators
    Z = compact_(Z,'zeros',eps);

    % check if reduction is required
    [n, nrOfGens] = size(generators(Z));

    if n*order < nrOfGens
        
        % reduce with the selected method
        if strcmp(method,'sum')
            Z = aux_reduceUnderApproxSum(Z,order);
        elseif strcmp(method,'scale')
            Z = aux_reduceUnderApproxScale(Z,order);
        elseif strcmp(method,'linProg')
            Z = aux_reduceUnderApproxLinProg(Z,order);
        elseif strcmp(method,'wetzlinger')
            Z = aux_reduceUnderApproxWetzlinger(Z,order);
        end
    else
        return;
    end

end


% Auxiliary functions -----------------------------------------------------

function Zred = aux_reduceUnderApproxLinProg(Z,order)
% reduce the zonotope order by computing an interval under-approximation of
% the zonotope spanned by the reduced generators using linear programming
   
    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N);

    % construct zonotope from the generators that are reduced
    Z1 = zonotope([zeros(length(c),1),Gred]);
    
    % construct zonotope from an interval enclosure
    Z2 = reduce(Z1,'pca',1);

    % compute largest interval that fits inside the zonotope by using
    % linear programming
    [A,b,Aeq,beq,lb,ub,f,ind] = aux_containmentConstraints(Z2,Z1);
    
    % solve linear program
    persistent options
    if isempty(options)
        options = optimoptions('linprog','Algorithm','interior-point', ...
                               'MaxIterations',10000,'display','off');
    end

    [s,~,exitflag] = linprog(f',A,b,Aeq,beq,lb,ub,options);
    
    if exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end
    
    s = s(ind);

    % construct the reduced zonotope
    Zred = zonotope([c,G,generators(Z2)*diag(s)]);

end

function Zred = aux_reduceUnderApproxScale(Z,order)
% An over-approximative reduced zonotope is computed first. This zonotope
% is then scaled using linear programming until it is fully contained in 
% the original zonotope 

    % over-approximative reduction of the zonotope
    Z_ = reduce(Z,'girard',order);
    
    % get conditions for linear program to scale the over-approximative 
    % zonotope until it is contained inside the original zonotope
    [A,b,Aeq,beq,lb,ub,f,ind] = aux_containmentConstraints(Z_,Z);
    
    % solve linear program
    persistent options
    if isempty(options)
        options = optimoptions('linprog','Algorithm','interior-point', ...
                               'MaxIterations',10000,'display','off');
    end             
                       
    [s,~,exitflag] = linprog(f',A,b,Aeq,beq,lb,ub,options);
    
    if exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end
    
    s = s(ind);
    
    % construct final zonotope
    Zred = zonotope([center(Z_),generators(Z_)*diag(s)]);
end

function Zred = aux_reduceUnderApproxSum(Z,order)
% sum up the generators that are reduced to obtain an inner-approximation

    % select generators to reduce
    n = dim(Z);
    N = floor(order*n - 1);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N);

    % under-approximate the generators that are reduced by one generator
    % corresponding to the sum of generators
    g = sum(Gred,2);

    % construct the reduced zonotope object
    Zred = zonotope([c,G,g]);
end

function Zred = aux_reduceUnderApproxWetzlinger(Z,order)
% reduction based on the Hausdorff distance betwenn a zonotope and its
% interval enclosure (see Theorem 3.2 in [2])

    % select generators to reduce
    n = dim(Z);
    N = floor(order*n-n);
    
    [c,G,Gred] = aux_selectSmallestGenerators(Z,N);

    % construct zonotope from the generators that are reduced
    Z1 = zonotope([zeros(length(c),1),Gred]);
    
    % state space transformation
    [S,~,~] = svd([-Gred,Gred]);
    Z1 = S' * Z1;
    
    % compute over-approximation of the Hausdorff distance between the
    % zonotope and its box enclsoure according to Theorem 3.2 in [2]
    G_ = abs(generators(Z1));
    
    for i = 1:size(G_,2)
        [~,ind] = max(G_(:,i));
        G(ind(1),i) = 0;
    end
    
    dist = 2*norm(sum(G_,2));
    
    % compute interval inner-approximation of the generators that are
    % reduced using the Minkowski difference for intervals
    int1 = interval(Z1);
    int2 = dist*interval(-ones(n,1),ones(n,1));
    
    int = minkDiff(int1,int2);
    
    % combine the interval inner-approximation with the unreduced
    % generators
    if ~isempty(int)
        Zred = zonotope(c,G) + S * zonotope(int);
    else
        Zred = zonotope(c,G); 
    end
end

function [c,G,Gred] = aux_selectSmallestGenerators(Z,N)
% select the generators that are reduced

    % obtain object properties
    c = center(Z);
    G_ = generators(Z);

    % sort according to generator length
    temp = sum(G_.^2,1);
    [~,ind] = sort(temp,'descend');
    
    % split into reduced and unreduced generators
    G = G_(:,ind(1:N));
    Gred = G_(:,ind(N+1:end));
end

function [A,b,Aeq,beq,lb,ub,f,ind] = aux_containmentConstraints(Zx,Zy)
% provides the contraints for linear programming for zonotope in zonotope 
% containment

    % decide depending on the complexity of the problem if an approximate
    % or the exact containment problem is solved
    comb = combinator(size(Zy.G,2),dim(Zy),'c');
    
    if size(comb,1) < 1000
        [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrPolytope(Zx,Zy);
    else
        [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrTedrake(Zx,Zy);
    end
end

function [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrPolytope(Zx,Zy)
% construct inequality constraints for zonotope X in zonotope Y containment

    % get halfspace representation of the zonotope
    poly = polytope(zonotope(Zy));
    C = poly.A;
    d = poly.b;

    % construct inequality constraints for all polytope halfspaces
    G = generators(Zx);
    c = center(Zx);
    
    A = zeros(size(C,1),size(G,2));
    b = zeros(size(C,1),1);
    
    for i = 1:size(C,1)
       A(i,:) = abs(C(i,:)*G);
       b(i) = d(i) - C(i,:)*c; 
    end
    
    % lower and upper bound
    ub = ones(size(G,2),1);
    lb = zeros(size(G,2),1);
    
    % objective function
    f = -sqrt(sum(G.^2,1));
    
    % equality constraints
    Aeq = [];
    beq = [];
    
    ind = 1:length(ub);

end

function [A,b,Aeq,beq,lb,ub,f,ind] = aux_contConstrTedrake(Zx,Zy)
% this function returns sufficient conditions for the zonotope X to be
% contained in zonotope Y according to Theorem 1 in [1]. The optimization 
% variable x is defined as follows:     
%
%   x = [s,T(:,1),...,T(:,nx),beta,Lambda(1,:),...,Lambda(qy,:)],
%
% where s \in [0,\inf] are the scaling factors and T,beta, and Lamda are
% auxiliary variables

    % get zonotope center and generator matrix
    Gx = generators(Zx); cx = center(Zx);
    Gy = generators(Zy); cy = center(Zy);

    % get required variables
    n = length(cx);
    
    nx = size(Gx,2);
    ny = size(Gy,2);
    
    qx = 2*nx; qy = 2*ny;
    
    Hx = [eye(nx);-eye(nx)];
    hx = [ones(nx,1);ones(nx,1)];
    
    Hy = [eye(ny);-eye(ny)];
    hy = [ones(ny,1);ones(ny,1)];

    % constraint Gx*diag(s) = Gy*T
    temp = repmat({Gy},[nx,1]);
    A1_ = blkdiag(temp{:});
    
    temp = num2cell(Gx,1);
    A2_ = blkdiag(temp{:});
    
    Aeq1 = [A2_,-A1_,zeros(size(A1_,1),ny + qx*qy)];
    beq1 = zeros(size(Aeq1,1),1);
    
    % constraint cy - cx = Gy*beta
    Aeq2 = [zeros(n,nx+nx*ny),Gy,zeros(n,qx*qy)];
    beq2 = cy - cx;
    
    % constraint Lambda*Hx = Hy*T
    temp = repmat({Hy},[nx,1]);
    A1_ = blkdiag(temp{:});
    
    temp = num2cell(Hx,1);
    A2_ = [];
    
    for i = 1:length(temp)
       temp_ = repmat({temp{i}'},[qy,1]);
       A2_ = [A2_;blkdiag(temp_{:})];
    end
    
    Aeq3 = [zeros(size(A1_,1),nx),-A1_,zeros(size(A1_,1),ny),A2_];
    beq3 = zeros(size(A1_,1),1);
    
    % constraint Lambda*hx < hy + Hy*beta
    temp = repmat({hx'},[qy,1]);
    A_ = blkdiag(temp{:});

    A1 = [zeros(size(A_,1),nx+nx*ny),-Hy,A_];
    b1 = hy;
    
    % constraint s > 0
    A2 = [-eye(nx),zeros(nx,size(A1,2)-nx)];
    b2 = zeros(nx,1);
    
    % constraint s < 1
    A3 = [eye(nx),zeros(nx,size(A1,2)-nx)];
    b3 = ones(nx,1);
    
    % constraint Lambda > 0
    A4 = [zeros(qx*qy,nx+nx*ny+ny),-eye(qx*qy)];
    b4 = zeros(qx*qy,1);
    
    % objective function
    f = zeros(size(A1,2),1);
    f(1:nx,1) = -sqrt(sum(Gx.^2,1))';
    
    % assemble overall constraint matrices
    Aeq = [Aeq1;Aeq2;Aeq3];
    beq = [beq1;beq2;beq3];
    
    A = [A1;A2;A3;A4];
    b = [b1;b2;b3;b4];
    
    lb = [];
    ub = [];
    
    ind = 1:nx;
end

% ------------------------------ END OF CODE ------------------------------
