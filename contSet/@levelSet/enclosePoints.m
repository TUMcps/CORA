function ls = enclosePoints(points,varargin)
% enclosePoints - enclose a point cloud with a level set
%
% Syntax:
%    ls = enclosePoints(points)
%    ls = enclosePoints(points,method)
%    ls = enclosePoints(points,method,order)
%
% Inputs:
%    points - matrix storing point cloud (dimension: [n,p] for p points)
%    method - method used for calculation
%               - 'multiple': level set with multiple constraints (default)
%               - 'single': level set consisting of a single constraint
%    order - max. polynomial order for polynomial constraints (default: 4)
%
% Outputs:
%    ls - levelSet object
%
% Example: 
%    p = -1 + 2*rand(2,10);
%
%    ls = levelSet.enclosePoints(p);
%    
%    figure; hold on; box on;
%    xlim([-1.5,1.5]); ylim([-1.5,1.5]);
%    plot(ls);
%    plot(p(1,:),p(2,:),'.k');
%
% References:
%    [1] A. A. Ahmadi and P. A. Parrilo. (2017). Sum of Squares 
%        Certificates for Stability of Planar, Homogeneous, and Switched 
%        Systems. IEEE Transactions on Automatic Control.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: levelSet

% Authors:       Niklas Kochdumper
% Written:       28-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    [method,order] = setDefaultValues({'multiple',4},varargin);

    % check input arguments
    inputArgsCheck({{points,'att','numeric','nonempty'};
                    {method,'str',{'single','multiple'}}; ...
                    {order,'att','numeric', ...
                            {'nonempty','scalar','integer','positive'}}});

    % compute enclosing level set with the selected method
    if strcmp(method,'single')
        ls = aux_enclosePointsSingle(points,order);
    elseif strcmp(method,'multiple')
        ls = aux_enclosePointsMultiple(points,order);
    end
end


% Auxiliary functions -----------------------------------------------------

function ls = aux_enclosePointsMultiple(points,order)
% compute an enclosing level set which consists of the intersection of
% multiple inequality constraints

    % center the point cloud at the origin
    n = size(points,1);
    c = mean(points,2);

    points = points - c;

    % construct all momomials for the given polynomial order
    comb = combinator(order+1,n)-1;
    comb = comb(sum(comb,2) <= order,:);

    % transform data points by the polynomial function
    data = zeros(size(points,2),size(comb,1));

    for i = 1:size(points,2)
        data(i,:) = prod(points(:,i).^(comb'),1);
    end

    % constraint a1 + a2*x1 + a3*x2 + a4*x1*x2 + ... < 0 that ensurese that
    % all points are located inside the level set
    A = data; b = zeros(size(A,1),1);

    % constraint dmax > a1 + a2*x1 + a3*x2 + a4*x1*x2 + ... that ensures
    % that variable dmax represents the maximum distance of a point to the
    % level set boundary
    A = [[zeros(size(A,1),1),A]; ...
         [-ones(size(A,1),1),-A]];
    
    b = [b;zeros(size(points,2),1)];

    % constraint a1 + a2 + a3 + ... = +/-1 (to avoid trivial solution
    % a1,a2,a3,... = 0)
    Aeq = ones(1,size(A,2)); Aeq(1) = 0;
    beq = 1;

    % minimize dmax (maximum distance of a point to the boundary) 
    f = zeros(size(A,2),1); f(1) = 1;

    % determinine optimal coefficients via linear programming
    problem.f = f; problem.Aineq = A; problem.bineq = b; 
    problem.Aeq = Aeq; problem.beq = beq;
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    problem.options = options;
    
    [x1,val1,exitflag1] = CORAlinprog(problem);
    
    % change the constraint constraint a1 + a2 + a3 + ... = 1 to 
    % a1 + a2 + a3 + ... = -1 and take the better solution of the two
    problem.beq = -problem.beq;
    [x2,val2,exitflag2] = CORAlinprog(problem);

    if exitflag1 < 0 && exitflag2 < 0
        throw(CORAerror('CORA:solverIssue'));
    end

    if exitflag2 < 0 || (exitflag1 > 0 && val1 < val2)
        x = x1;
    else
        x = x2;
    end

    coeff = x(2:end);

    % construct level set object
    x = sym('x',[n,1]);
    eq = 0;

    for i = 1:length(coeff)
        eq = eq + coeff(i)*prod((x-c).^(comb(i,:)'));
    end

    ls = levelSet(eq,x,'<=');

    % intersect with enclosing interval
    lsInt = levelSet(polytope(interval.enclosePoints(points) + c));

    ls = ls & lsInt;
end

function ls = aux_enclosePointsSingle(points,order)
% compute an enclosing level set which consists of a single inequality
% constraint

    % initialization
    order = max(2,order);

    if order > 1 && mod(order,2) ~= 0
        order = order - 1;
    end

    % center the point cloud at the origin
    n = size(points,1);
    c = mean(points,2);

    points = points - c;

    % construct all momomials for the given polynomial order
    comb = combinator(order+1,n)-1;
    comb = comb(sum(comb,2) <= order,:);

    % remove non-square monomials of highest order to ensure highest order
    % monomials are positive definite (required for the level set to be 
    % bounded, see paragraph after Proposition IV.1 in [1])
    comb = comb(sum(comb,2) <= order | max(mod(comb,2),[],2) == 0,:);

    % transform data points by the polynomial function
    data = zeros(size(points,2),size(comb,1));

    for i = 1:size(points,2)
        data(i,:) = prod(points(:,i).^(comb'),1);
    end

    % constraint a1 + a2*x1 + a3*x2 + a4*x1*x2 + ... < 0 that ensurese that
    % all points are located inside the level set
    A = data; b = zeros(size(A,1),1);

    % constraint dmax > a1 + a2*x1 + a3*x2 + a4*x1*x2 + ... that ensures
    % that variable dmax represents the maximum distance of a point to the
    % level set boundary
    A = [[zeros(size(A,1),1),A]; ...
         [-ones(size(A,1),1),-A]];
    
    b = [b;zeros(size(data,1),1)];

    % constraint highest-degree polynomial positive definite (required for
    % level set to be bounded, see paragraph after Proposition IV.1 in [1])
    ind = find(sum(comb,2) == order);

    A_ = zeros(length(ind),size(A,2));
    A_(:,ind+1) = -eye(length(ind));
    b_ = -0.1*ones(length(ind),1);

    A = [A;A_]; b = [b;b_];

    % minimize dmax (maximum distance of a point to the boundary) 
    f = zeros(size(A,2),1); f(1) = 1;

    % determinine optimal coefficients via linear programming
    problem.f = f; problem.Aineq = A; problem.bineq = b; 
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    problem.options = options;
    
    [x,~,exitflag] = CORAlinprog(problem);
    
    if exitflag < 0
        throw(CORAerror('CORA:solverIssue'));
    end

    coeff = x(2:end);

    % construct level set object
    x = sym('x',[n,1]);
    eq = 0;

    for i = 1:length(coeff)
        eq = eq + coeff(i)*prod((x-c).^(comb(i,:)'));
    end

    ls = levelSet(eq,x,'<=');
end

% ------------------------------ END OF CODE ------------------------------
