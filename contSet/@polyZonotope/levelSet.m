function ls = levelSet(pZ,varargin)
% levelSet - computes a level set which encloses the polynomial zonotope
%
% Syntax:
%    ls = levelSet(pZ)
%    ls = levelSet(pZ,order)
%
% Inputs:
%    pZ - polyZonotope object
%    order - maximum order of polynomials for the level set (default: 2)
%
% Outputs:
%    ls - levelSet object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1 1;0 2 1 2],[0;0],[1 0 3 1;0 1 0 2]);
%    ls = levelSet(pZ);
%
%    figure; hold on;
%    xlim([-5,5]); ylim([-5,5]);
%    plot(ls);
%    plot(pZ,[1,2],'LineWidth',2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Niklas Kochdumper
% Written:       30-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    order = setDefaultValues({2},varargin);

    % check input arguments
    inputArgsCheck({{order,'att','numeric', ...
                            {'nonempty','scalar','integer','positive'}}});

    % initialization
    n = dim(pZ);
    ls = fullspace(n);

    % determine the n most important factors
    [ind,pZ_] = aux_mostImportantFactors(pZ,n);

    f = fhandle(pZ_); fac = min(n,length(ind)); ind = ind(1:fac);

    % generate uniformly distributed points on the inside and on the
    % boundary of the polynomial zonotope
    Iinside = interval(-ones(fac,1),ones(fac,1));
    Iboundary = interval(-ones(max(1,fac-1),1),ones(max(1,fac-1),1));

    aInside = aux_uniformDistributedPoints(Iinside,1000);
    aBoundary = aux_uniformDistributedPoints(Iboundary,1000);

    pInside = zeros(n,size(aInside,2));

    for i = 1:size(aInside,2)
        pInside(:,i) = f(aInside(:,i));
    end

    % loop over all (n-1)-dimensional faces of the hypercube spanned by the
    % most important factors
    tmp = cellfun(@(x) x*[1,1],num2cell(1:fac),'UniformOutput',false);
    F = [repmat([1,-1],[1,fac]);[tmp{:}]];

    for i = 1:size(F,2)
    
        % get points from the corresponding boundary
        if fac == 1
            a = aBoundary;
        else
            a = F(1,i)*ones(fac,size(aBoundary,2));
            a(ind(ind ~= F(2,i)),:) = aBoundary;
        end

        p = zeros(n,size(aBoundary,2));

        for j = 1:size(aBoundary,2)
            p(:,j) = f(a(:,j));
        end

        % construct template level set by fitting the boundary points
        val = inf;

        for o = 1:order

            % construct template level set by fitting the boundary points
            lsTmp = aux_templateLevelSet(p,pInside,o);

            % shift level set to guarantee it contains the poly. zonotope
            if isempty(symvar(lsTmp.eq))
                I = interval(eval(lsTmp.eq));
            else
                tayEq = lsTmp.funHan(taylm(pZ));
                I = interval(tayEq,'bnbAdv');
            end

            ls_ = levelSet(lsTmp.eq - supremum(I),lsTmp.vars,'<=');

            % compute accuracy of the level set
            val_ = sum(abs(ls_.funHan(p)));

            if val_ < val
                val = val_; lsNew = ls_;
            end
        end

        % intersect with previous level set
        ls = ls & lsNew;
    end

    % intersect with enclosing interval to ensure that the level set is
    % bounded
    ls = ls & levelSet(polytope(interval(pZ)));
end


% Auxiliary functions -----------------------------------------------------

function [ind,pZ_] = aux_mostImportantFactors(pZ,n)
% construct a modifed polynomial zonotope that only consists of the n most
% important factors

    p = size(pZ.E,1);

    % sort factors according to their importance (= length of the
    % corresponding generators)
    isEvenColumn = all(mod(pZ.E, 2) == 0, 1);

    G = pZ.G; G(:,isEvenColumn) = 0.5*G(:,isEvenColumn);
    lenG = sum(G.^2,1); len = zeros(p,1);

    for i = 1:p
        len(i) = sum(lenG(pZ.E(i,:) > 0));
    end

    [~,ind] = sort(len,'descend');

    % different cases depending on the number of factors
    if p == n

        pZ_ = polyZonotope(pZ.c,pZ.G,[],pZ.E(ind,:));

    elseif p <= n

        lenGI = sum(pZ.GI.^2,1);
        [~,ind_] = sort(lenGI,'descend');
        ind_ = ind_(1:min(length(ind_),n-p));

        pZ_ = polyZonotope(pZ.c,[pZ.G,pZ.GI(:,ind_)],[], ...
                                   blkdiag(pZ.E(ind,:),eye(length(ind_))));
        ind = 1:size(pZ_.E,1);

    else
        pZ_ = polyZonotope(pZ.c,pZ.G,[],pZ.E(ind(1:n),:));
        ind = 1:n;
    end
end

function p = aux_uniformDistributedPoints(I,N)
% generate uniformly distributed points inside an interval

    % get properties from the interval
    c = center(I); lb = I.infimum; ub = I.supremum;

    % determine how many points to generate for each dimension
    n = dim(I); num = ones(n,1); cnt = 1;

    while prod(num) < N
        num(mod(cnt,n)+1) = num(mod(cnt,n)+1) + 1;
        cnt = cnt + 1;
    end

    % generate points for each dimension
    pDim = cell(length(num),1);

    for i = 1:length(num)
        if num(i) == 1
            pDim{i} = c(i);
        else
            pDim{i} = linspace(lb(i),ub(i),num(i));
        end
    end

    % combine points for different dimension via Cartesian product
    p = pDim{1};

    for i = 2:length(num)
        
        p_ = [];

        for j = 1:num(i)
            p_ = [p_,[p;ones(1,size(p,2))*pDim{i}(j)]];
        end

        p = p_;
    end
end

function ls = aux_templateLevelSet(pBoundary,pInside,order)
% construct a template level set by fitting the boundary points

    % combine boundary points and points from the inside
    points = [pBoundary,pInside];
    ind = 1:size(pBoundary,2);

    % center the point cloud at the origin
    n = size(points,1);
    c = mean(pBoundary,2);

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
         [-ones(length(ind),1),-A(ind,:)]];
    
    b = [b;zeros(length(ind),1)];

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
end

% ------------------------------ END OF CODE ------------------------------
