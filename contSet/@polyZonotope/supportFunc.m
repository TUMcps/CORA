function val = supportFunc(pZ,dir,varargin)
% supportFunc - Calculates the upper or lower bound of a polynomial
%    zonotope along a given direction
%
% Syntax:  
%    val = supportFunc(pZ,dir)
%    val = supportFunc(pZ,dir,type)
%    val = supportFunc(pZ,dir,type,method)
%    val = supportFunc(pZ,dir,type,'bnb',maxOrder)
%    val = supportFunc(pZ,dir,type,'bnbAdv',maxOrder)
%    val = supportFunc(pZ,dir,type,'globOpt',maxOrder,tol)
%    val = supportFunc(pZ,dir,type,'split',splits)
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
%    tol - tolerance for the verified global optimization method
%    split - number of splits that are performed to calculate the bounds
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
% See also: interval, conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      29-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % pre-processing
    [res,vars] = pre_supportFunc('polyZonotope',pZ,dir,varargin{:});
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of var
        val = vars{1}; return
    else
        pZ = vars{1}; dir = vars{2}; type = vars{3};
        method = vars{4}; tol = vars{6};
        if strcmp(method,'split')
            splits = vars{5};
        elseif any(strcmp(method,{'bnb','bnbAdv','globOpt'}))
            maxOrder = vars{5};
        end
    end
    

    % different methods to evaluate the support function
    if strcmp(method,'interval')
       
        % project the polynomial zonotope onto the direction
        pZ_ = dir' * pZ;
        
        % compute support function based on interval enclosure
        val = interval(zonotope(pZ_));
        
        if strcmp(type,'lower')
            val = infimum(val); 
        elseif strcmp(type,'upper')
            val = supremum(val);
        elseif strcmp(type,'range')
            % 'val' is already desired result
        end
        
    elseif strcmp(method,'bernstein')
        
        val = supportFuncBernstein(pZ,dir,type);
        
    elseif strcmp(method,'split')
        
        val = supportFuncSplit(pZ,dir,type,splits);

    % determine bounds with taylor models and "branch and bound" or
    % advanced "branch and bound"
    elseif strcmp(method,'bnb') || strcmp(method,'bnbAdv')

        val = supportFuncBnB(pZ,dir,type,method,maxOrder);

    % determine bounds with global verified optimization
    elseif strcmp(method,'globOpt')

        val = supportFuncGlobOpt(pZ,dir,type,maxOrder,tol);

    % determine bounds using quadratic programming
    elseif strcmp(method,'quadProg')

        val = supportFuncQuadProg(pZ,dir,type);

    end
end


% Auxlilary functions -----------------------------------------------------

function val = supportFuncBernstein(pZ,dir,type)
% compute the support function using Bernstein polynomials

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;

    % initialization
    p = size(pZ_.expMat,1);    
    dom = interval(-ones(p,1),ones(p,1));  

    % dependent generators: convert to bernstein polynomial
    B = poly2bernstein(pZ_.G,pZ_.expMat,dom);

    int1 = interval(min(min(B)),max(max(B)));

    % independent generators: enclose zonotope with interval
    int2 = interval(zonotope([pZ_.c,pZ_.Grest]));

    val = int1 + int2;

    if strcmp(type,'lower')
       val = infimum(val); 
    elseif strcmp(type,'upper')
       val = supremum(val);
    end
end

function val = supportFuncSplit(pZ,dir,type,splits)
% compute support function by recursively splitting the set

    % handle different types
    if strcmp(type,'lower')
       val = -supportFuncSplit(pZ,-dir,'upper',splits);
       return;
    elseif strcmp(type,'range')
       up = supportFuncSplit(pZ,dir,'upper',splits);
       low = -supportFuncSplit(pZ,-dir,'upper',splits);
       val = interval(low,up);
       return;
    end

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;

    % split the polynomial zonotope multiple times to obtain a better 
    % over-approximation of the real shape
    pZsplit{1} = pZ_;
    val_min = -inf; val_max = -inf;

    for i = 1:splits
        
        qZnew = [];
        
        for j = 1:length(pZsplit)
            
            res = splitLongestGen(pZsplit{j});
            
            for k = 1:length(res)
                
                % compute support function for enclosing zonotope
                [val,~,alpha] = supportFunc(zonotope(res{k}),1);
                
                % update upper and lower bound
                if val > val_max
                   
                    % update upper bound
                    val_max = val;
                    
                    % update lower bound by determining most critical point
                    ind1 = find(sum(res{k}.expMat,1) == 1);
                    ind2 = find(sum(res{k}.expMat(:,ind1),2) == 1);
                    
                    beta = alpha(size(res{k}.expMat,2)+1:end);
                    alpha_ = zeros(size(res{k}.expMat,1),1);
                    alpha_(ind2) = alpha(ind1);
                    
                    tmp = res{k}.c + sum(res{k}.G .*  ...
                                            prod(alpha_.^res{k}.expMat,1)); 
                    
                    if ~isempty(beta)
                        tmp = tmp + res{k}.Grest*beta;
                    end
                    
                    if tmp > val_min
                       val_min = tmp; 
                    end
                end
                
                % add new set to queue (if not smaller than lower bound)
                if val > val_min
                   qZnew{end+1} = res{k}; 
                end
            end
        end
        
        pZsplit = qZnew;
    end
    
    val = val_max;
end

function val = supportFuncBnB(pZ,dir,type,method,maxOrder)
% compute the support function using branch and bound methods

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;
    
    % over-approximate the independent generators
    valInd = interval(zonotope([pZ_.c,pZ_.Grest]));

    % create taylor models
    n = size(pZ_.expMat,1);
    temp = ones(n,1);
    T = taylm(interval(-temp,temp),maxOrder,[],method);

    % calculate bounds of the polynomial part (= dependent
    % generators)
    valDep = interval(polyPart(T,pZ_.G,pZ_.expMat));

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

function val = supportFuncGlobOpt(pZ,dir,type,maxOrder,tol)
% compute the support function using verified glboal optimization

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;

    % over-approximate the independent generators
    valInd = interval(zonotope([pZ_.c,pZ_.Grest]));

    % remove zero entries from the generator matrix
    ind = find(pZ_.G == 0);
    G = pZ_.G;
    expMat = pZ_.expMat;
    G(:,ind) = [];
    expMat(:,ind) = [];

    % remove zero rows from the exponent matrix
    expMat(sum(expMat,2) == 0,:) = [];

    % function to opimize
    func = @(x) polyPart(x,G,expMat);

    % domain for the optimization variables
    n = size(expMat,1);
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

function val = supportFuncQuadProg(pZ,dir,type)
% compute the support function using quadratic programming

    % consider different types of bounds
    if strcmp(type,'range')
        l = supportFunc(pZ,dir,'lower');
        u = supportFunc(pZ,dir,'upper');
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
    indLin = find(sum(pZ_.expMat,1) <= 1);
    p = length(pZ_.id); f = zeros(p,1);

    for i = 1:length(indLin)
        ind = find(pZ_.expMat(:,indLin(i)) == 1);
        f(ind(1)) = pZ_.G(:,indLin(i));
    end

    % extract quadratic part x'*H*x for states
    indQuad = find(sum(pZ_.expMat,1) == 2);
    H = zeros(p);

    for i = 1:length(indQuad)
       if max(pZ_.expMat(:,indQuad(i))) == 2
           ind = find(pZ_.expMat(:,indQuad(i)) == 2);
           H(ind(1),ind(1)) = pZ_.G(:,indQuad(i));
       else
           ind = find(pZ_.expMat(:,indQuad(i)) == 1);
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
        val = supportFunc(pZ,dir,type,'interval');
        return;
    end
    if ~isempty(ind_)
        temp = zeros(size(d)); temp(ind_) = d(ind_);
        N = N + V*diag(temp)*V';
    end

    % enclose remaining part with additional factors
    ind = setdiff(1:size(pZ_.G,2),[indQuad,indLin]);
    G = pZ_.G(:,ind); expMat = pZ_.expMat(:,ind);

    for i = 1:size(N,1)
        for j = i:size(N,2)
            if i == j && N(i,j) ~= 0
                G = [G, N(i,j)];
                temp = zeros(p,1); temp(i) = 2;
                expMat = [expMat,temp];
            elseif N(i,j) ~= 0 || N(j,i) ~= 0
                G = [G, N(i,j) + N(j,i)];
                temp = zeros(p,1); temp(i) = 1; temp(j) = 1;
                expMat = [expMat,temp];
            end
        end
    end

    % extend matrices of quadratic program
    if ~isempty(expMat)
        Z = zonotope(polyZonotope(c,G,[],expMat));

        H = blkdiag(H,zeros(size(Z.Z,2)-1));
        f = [f; generators(Z)']; c = center(Z);
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
    if ~isempty(pZ.Grest)
       val = val + supportFunc(zonotope(zeros(dim(pZ),1), ...
                               pZ.Grest),dir,type); 
    end
end

function val = polyPart(x,G,expMat)

    for i = 1:length(G)
       temp = x(1)^expMat(1,i);
       for j = 2:size(expMat,1)
           temp = temp * x(j)^expMat(j,i);
       end
       if i == 1
           val = G(i) * temp;
       else
           val = val + G(i) * temp;
       end
    end
end


%------------- END OF CODE --------------