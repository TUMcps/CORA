function val = supportFunc_(pZ,dir,type,method,maxOrderOrSplits,tol,varargin)
% supportFunc_ - Calculates the upper or lower bound of a polynomial
%    zonotope along a given direction
%
% Syntax:
%    val = supportFunc_(pZ,dir)
%    val = supportFunc_(pZ,dir,type)
%    val = supportFunc_(pZ,dir,type,method)
%    val = supportFunc_(pZ,dir,type,'bnb',maxOrder)
%    val = supportFunc_(pZ,dir,type,'bnbAdv',maxOrder)
%    val = supportFunc_(pZ,dir,type,'globOpt',maxOrder,tol)
%    val = supportFunc_(pZ,dir,type,'split',splits)
%
% Inputs:
%    pZ - polyZonotope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    method - method that is used to calculate the bounds for the dependent
%             part of the polynomial zonotope
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'bnb': taylor models with "branch and bound" algorithm
%              'bnbAdv': taylor models with advandced bnb-algorithm
%              'globOpt': verified global optimization 
%              'bernstein': conversion to a bernstein polynomial
%              'quadProg': quadratic programming
%    maxOrder - maximum polynomial order of the taylor model
%    splits - number of splits that are performed to calculate the bounds
%    tol - tolerance for the verified global optimization method
%
% Outputs:
%    val - interval object specifying the upper and lower bound along the
%          direction
%
% Examples:
%    pZ = polyZonotope([0;0],[2 -1 2;0 -2 -3],[],[1 0 1;0 1 3]);
%    dir = [1;2];
%    
%    b1 = supportFunc(pZ,dir,'range','interval');
%    ub1 = conHyperplane(dir,supremum(b1));
%    lb1 = conHyperplane(dir,infimum(b1));
%
%    b2 = supportFunc(pZ,dir,'range','globOpt');
%    ub2 = conHyperplane(dir,supremum(b2));
%    lb2 = conHyperplane(dir,infimum(b2));
%
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(ub1,[1,2],'b');
%    plot(lb1,[1,2],'b');
%    plot(ub2,[1,2],'g');
%    plot(lb2,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, interval, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       29-July-2018
% Last update:   17-October-2022 (NK, improve 'split' method)
%                06-December-2022 (TL, fix: 'split' considers splitted sets)
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

    if strcmp(method,'split')
        splits = maxOrderOrSplits;
    elseif any(strcmp(method,{'bnb','bnbAdv','globOpt'}))
        maxOrder = maxOrderOrSplits;
    end

    % different methods to evaluate the support function
    if strcmp(method,'interval')
       
        % project the polynomial zonotope onto the direction
        pZ_ = dir' * pZ;
        
        % compute support function based on interval enclosure
        val = interval(zonotope(pZ_));
        
        if strcmp(type,'lower')
            if representsa_(val,'emptySet',eps)
                val = Inf;
            else
                val = infimum(val);
            end            
        elseif strcmp(type,'upper')
            if representsa_(val,'emptySet',eps)
                val = -Inf;
            else
                val = supremum(val);
            end
        elseif strcmp(type,'range')
            if representsa_(val,'emptySet',eps)
                val = [-Inf,Inf];
            end
            % otherwise 'val' is already desired result
        end
        
    elseif strcmp(method,'bernstein')
        
        val = aux_supportFuncBernstein(pZ,dir,type);
        
    elseif strcmp(method,'split')
        
        val = supportFuncSplit(pZ,dir,type,splits);

    % determine bounds with taylor models and "branch and bound" or
    % advanced "branch and bound"
    elseif strcmp(method,'bnb') || strcmp(method,'bnbAdv')

        val = aux_supportFuncBnB(pZ,dir,type,method,maxOrder);

    % determine bounds with global verified optimization
    elseif strcmp(method,'globOpt')

        val = aux_supportFuncGlobOpt(pZ,dir,type,maxOrder,tol);

    % determine bounds using quadratic programming
    elseif strcmp(method,'quadProg')

        val = aux_supportFuncQuadProg(pZ,dir,type);

    end
end


% Auxiliary functions -----------------------------------------------------

function val = aux_supportFuncBernstein(pZ,dir,type)
% compute the support function using Bernstein polynomials

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;

    % initialization
    p = size(pZ_.E,1);    
    dom = interval(-ones(p,1),ones(p,1));  

    % dependent generators: convert to bernstein polynomial
    B = poly2bernstein(pZ_.G,pZ_.E,dom);

    int1 = interval(min(min(B)),max(max(B)));

    % independent generators: enclose zonotope with interval
    int2 = interval(zonotope([pZ_.c,pZ_.GI]));

    val = int1 + int2;

    if strcmp(type,'lower')
       val = infimum(val); 
    elseif strcmp(type,'upper')
       val = supremum(val);
    end
end

function val = aux_supportFuncBnB(pZ,dir,type,method,maxOrder)
% compute the support function using branch and bound methods

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;
    
    % over-approximate the independent generators
    valInd = interval(zonotope([pZ_.c,pZ_.GI]));

    % create taylor models
    n = size(pZ_.E,1);
    temp = ones(n,1);
    T = taylm(interval(-temp,temp),maxOrder,[],method);

    % calculate bounds of the polynomial part (= dependent
    % generators)
    valDep = interval(aux_polyPart(T,pZ_.G,pZ_.E));

    % add the two parts from the dependent and independent
    % generators
    val = valDep + valInd;

    % extract the desired bound
    if strcmp(type,'lower')
       val = infimum(val); 
    elseif strcmp(type,'upper')
       val = supremum(val);
    end
end

function val = aux_supportFuncGlobOpt(pZ,dir,type,maxOrder,tol)
% compute the support function using verified glboal optimization

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;

    % over-approximate the independent generators
    valInd = interval(zonotope([pZ_.c,pZ_.GI]));

    % remove zero entries from the generator matrix
    ind = find(pZ_.G == 0);
    G = pZ_.G;
    E = pZ_.E;
    G(:,ind) = [];
    E(:,ind) = [];

    % remove zero rows from the exponent matrix
    E(sum(E,2) == 0,:) = [];

    % function to opimize
    func = @(x) aux_polyPart(x,G,E);

    % domain for the optimization variables
    n = size(E,1);
    temp = ones(n,1);
    dom = interval(-temp,temp);

    % calculate the bounds of the depenedent part
    if strcmp(type,'lower')
        minDep = globVerMinimization(func,dom,tol,maxOrder);
        val = minDep + infimum(valInd);
    elseif strcmp(type,'upper')
        maxDep = globVerMinimization(@(x) -func(x),dom,tol,maxOrder);
        val = -maxDep + supremum(valInd);
    else
        valDep = globVerBounds(func,dom,tol,maxOrder);
        val = valDep + valInd;
    end
end

function val = aux_supportFuncQuadProg(pZ,dir,type)
% compute the support function using quadratic programming

    % consider different types of bounds
    if strcmp(type,'range')
        l = supportFunc_(pZ,dir,'lower','interval',8,1e-3);
        u = supportFunc_(pZ,dir,'upper','interval',8,1e-3);
        val = interval(l,u);
        return;
    end
    
    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;
    
    c = pZ_.c; pZ_.c = 0;
    if strcmp(type,'upper')
        pZ_ = (-1)*pZ_;
    end

    % extract linear part f*x 
    indLin = find(sum(pZ_.E,1) <= 1);
    p = length(pZ_.id); f = zeros(p,1);

    for i = 1:length(indLin)
        ind = find(pZ_.E(:,indLin(i)) == 1);
        f(ind(1)) = pZ_.G(:,indLin(i));
    end

    % extract quadratic part x'*H*x for states
    indQuad = find(sum(pZ_.E,1) == 2);
    H = zeros(p);

    for i = 1:length(indQuad)
       if max(pZ_.E(:,indQuad(i))) == 2
           ind = find(pZ_.E(:,indQuad(i)) == 2);
           H(ind(1),ind(1)) = pZ_.G(:,indQuad(i));
       else
           ind = find(pZ_.E(:,indQuad(i)) == 1);
           H(ind(1),ind(2)) = 0.5 * pZ_.G(:,indQuad(i));
           H(ind(2),ind(1)) = 0.5 * pZ_.G(:,indQuad(i));
       end
    end

    % split matirx H for the quadratic part into a positive definite
    % matrix and a remainder matrix
    M = 0.5*(H + H'); N = 0.5*(H - H');     % get symmetric matrix
    [V,D] = eig(M);  d = diag(D);           % compute eigenvalues   
    ind = find(d >= 0);                     % get positive eigenvalues
    ind_ = setdiff(1:length(d),ind);        % get negative eigenvalues

    if ~isempty(ind)
        temp = zeros(size(d)); temp(ind) = d(ind);
        H = V*diag(temp)*V'; 
    else
        val = supportFunc_(pZ,dir,type,'interval',8,1e-3);
        return;
    end
    if ~isempty(ind_)
        temp = zeros(size(d)); temp(ind_) = d(ind_);
        N = N + V*diag(temp)*V';
    end

    % enclose remaining part with additional factors
    ind = setdiff(1:size(pZ_.G,2),[indQuad,indLin]);
    G = pZ_.G(:,ind); E = pZ_.E(:,ind);

    for i = 1:size(N,1)
        for j = i:size(N,2)
            if i == j && N(i,j) ~= 0
                G = [G, N(i,j)];
                temp = zeros(p,1); temp(i) = 2;
                E = [E,temp];
            elseif N(i,j) ~= 0 || N(j,i) ~= 0
                G = [G, N(i,j) + N(j,i)];
                temp = zeros(p,1); temp(i) = 1; temp(j) = 1;
                E = [E,temp];
            end
        end
    end

    % extend matrices of quadratic program
    if ~isempty(E)
        Z = zonotope(polyZonotope(c,G,[],E));

        H = blkdiag(H,zeros(size(Z.G,2)));
        f = [f; Z.G']; c = Z.c;
    end

    % solve quadratic program
    options = optimoptions('quadprog','Display','off');
    w = warning();
    warning('off');

    [~,val] = quadprog(2*H,f,[],[],[],[],-ones(length(f),1), ...
                       ones(length(f),1),[],options);

    warning(w);

    val = val + c;

    if strcmp(type,'upper')
       val = -val; 
    end

    % compute bound for independent generators
    if ~isempty(pZ.GI)
       val = val + supportFunc_( ...
           zonotope(zeros(dim(pZ),1), pZ.GI), dir, type); 
    end
end

function val = aux_polyPart(x,G,E)

    for i = 1:length(G)
       temp = x(1)^E(1,i);
       for j = 2:size(E,1)
           temp = temp * x(j)^E(j,i);
       end
       if i == 1
           val = G(i) * temp;
       else
           val = val + G(i) * temp;
       end
    end
end


% ------------------------------ END OF CODE ------------------------------
