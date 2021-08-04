function res = in(pZ,obj,varargin)
% in - determines if obj is contained in the polynomial zonotope pZ
%
% Syntax:  
%    res = in(pZ,obj)
%    res = in(pZ,obj,type)
%
% Inputs:
%    pZ - polyZonotope object
%    obj - contSet object or point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - boolean whether obj is contained in pZ, or not
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);
%
%    p1 = [1;1];
%    p2 = [-1;3];
%    obj1 = polyZonotope([0;0],0.3*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    obj2 = polyZonotope([0;0],0.4*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
% 
%    in(pZ,p1,'approx')
%    in(pZ,p2,'approx')
%    in(pZ,obj1,'approx')
%    in(pZ,obj2,'approx')
% 
%    figure; hold on;
%    plot(pZ,[1,2],'b');
%    plot(p1(1),p1(2),'.g','MarkerSize',20);
%    plot(p2(1),p2(2),'.r','MarkerSize',20);
%    plot(obj1,[1,2],'g');
%    plot(obj2,[1,2],'r');
%
% References: 
%   [1] O. Mullier and et al. "General Inner Approximation of Vector-valued 
%       Functions", Reliable Computing, 2013
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/in, conZonotope/in
%
% Author:       Niklas Kochdumper
% Written:      13-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = false;

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && strcmp(varargin{1},'approx')
        type = 'approx';
    end
    
    % check user inputs 
    if strcmp(type,'exact')
        error(errNoExactAlg(pZ,obj));
    end
    
    % initialize variables
    fHan = @(x) funcPoly(x,pZ.c,pZ.G,pZ.Grest,pZ.expMat);
    jacHan = funHanJacobian(pZ.G,pZ.Grest,pZ.expMat);
        
    temp = ones(length(pZ.id) + size(pZ.Grest,2),1);
    X = interval(-temp,temp);
        
    % point in polynomial zonotope containment
    if isnumeric(obj)
        
        % try to prove that the point is inside the polynomial zonotope
        % using the approach from [1]
        Y = interval(obj);
        
        x = getFactorDomain(fHan,Y,X);
        X_ = interval(x - 1e-10* ones(size(x)), x + 1e-10* ones(size(x)));
        
        res = containsInterval(fHan,jacHan,X,X_,Y);    
        
    % set in zonotope containment
    elseif isa(obj,'interval') || isa(obj,'zonotope') || ...
           isa(obj,'mptPolytope') || isa(obj,'zonoBundle') || ...
           isa(obj,'conZonotope') || isa(obj,'taylm') || ...
           isa(obj,'polyZonotope') || isa(obj,'conPolyZono')
        
       % convert set to polynomial zonotope
       if ~isa(obj,'conPolyZono')
           obj = polyZonotope(obj);
       end
       
       % try to prove that set containment does not hold
       temp = disproveContainment(pZ,obj);
       
       if temp == 1
          return; 
       end
       
       % try to prove that set containment holds
       res = proveContainment(obj,fHan,jacHan,X);

    else
        error(noops(pZ,obj));
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = disproveContainment(pZ1,pZ2)
% use contraction to prove that pZ2 is not a subset of pZ1

    % construct polynomial constraint for the intersection
    c = pZ1.c - pZ2.c;
    G = [pZ1.G -pZ2.G];
    Grest = [pZ1.Grest,-pZ2.Grest];
    expMat = blkdiag(pZ1.expMat,pZ2.expMat);
        
    % contract the factor domain \alpha_k in [-1,1] based on the
    % polynomial constraint
    n = size(expMat,1) + size(Grest,2);
    dom = interval(-ones(n,1),ones(n,1));

    int = contractPoly(c,G,Grest,expMat,dom,'all',3);

    % pZ2 is proven to be not contained in pZ1 -> res = 1
    temp = int(size(expMat,1)+1:end);
    res = 0;
    
    if isempty(temp) || any(supremum(temp) < 1) || any(infimum(temp) > -1)
       res = 1; 
    end
end

function res = proveContainment(obj,fHan,jacHan,X)
% try to prove that the set is contained in the polynomial zonotope by
% recursiveley splitting the set and using the method in [1] to prove that
% an interval enclosure of the splitted sets is contained in the polynomial
% zonotope

    % maximum number of recursive splits
    splits = 8;

    % list storing the splitted sets that still need to be proven to be
    % contained in the polynomial zonotope
    list = {obj};
    cnt = 1;
    
    % loop over all splits
    while cnt <= splits
       
        list_ = {};
        
        % loop over all sets in the list
        for j = 1:length(list)
           
            % try to prove that an interval enclosure of the set is
            % contained in the polynomial zonotope using the method in [1]
            Y = interval(list{j});
        
            x = getFactorDomain(fHan,Y,X);
            X_ = x + pinv(jacHan(x))*(Y - center(Y));
        
            temp = containsInterval(fHan,jacHan,X,X_,Y); 
            
            % split the set if it is not contained
            if ~temp
               sets = splitLongestGen_(list{j});
               
               list_{end+1} = sets{1};
               list_{end+1} = sets{2};
            end
        end
        
        % increase number of splits if the number of sets decreases
        if cnt == splits && ~isempty(list_) && length(list_) <= length(list)
            splits = splits + 1;
        end
        
        % update list
        list = list_;
        cnt = cnt + 1;
    end
    
    % check if all splitted sets are proven to be contained
    res = isempty(list_);
end

function f = funcPoly(x,c,G,Grest,expMat)
% evaluate the mapping function of a polynomial zonotope

    % initialization
    f = c;
    n = size(expMat,1);
    x1 = x(1:n); x2 = x(n+1:end);
    
    % loop over all dependent generators
    for i = 1:size(G,2)
        f = f + G(:,i)*prod(x1.^expMat(:,i),1); 
    end

    % add independent generators
    if ~isempty(Grest)
        f = f + Grest*x2;
    end
end

function J = jacobianPoly(x,G,Grest,expMat,n,p)
% compute the jacobian matrix for the mapping function of a polynomial
% zonotope

    % initialization
    J = 0 * repmat(x(1),[n,length(expMat)]);
    x_ = x(1:p);
    
    % loop over all variables
    for i = 1:length(expMat)
       for j = 1:size(expMat{i},2)
          J(:,i) = J(:,i) + G{i}(:,j)*prod(x_.^expMat{i}(:,j),1); 
       end
    end
    
    % consider independent generators
    J = [J,Grest];
end

function jacHan = funHanJacobian(G,Grest,expMat)
% compute a function handle for the jacobian matrix for the mapping 
% function of a polynomial zonotope

    % compute exponent matrix differentiazed for each variable
    Elist = cell(size(expMat,1),1);
    Glist = cell(size(expMat,1),1);

    for i = 1:length(Elist)
       ind = find(expMat(i,:) > 0);
       temp = expMat(:,ind);
       temp(i,:) = temp(i,:) - 1;
       Elist{i} = temp;
       Glist{i} = G(:,ind) * diag(expMat(i,ind));
    end

    % function handle for jacobian matrix
    jacHan=@(x)jacobianPoly(x,Glist,Grest,Elist,size(G,1),size(expMat,1));
end

function res = containsInterval(f,df,X,X_,Y_,varargin)
% implementation of the function Inner() according to Alg. 3 in [2] which
% checks if an interval is part of the image of a nonlinear function

    % parameter
    tau = 1.01;
    mu = 0.9;
    
    % dermine suitable dimensions of the domain X
    J = df(X_);
    x_ = center(X_);
    n = size(J,1);
    
    ind = getSuitableSubmatrix(J,n);
    
    % check if Jacobian matrix has full rank
    if ismember(0,ind)  
        res = 0;
        return;
    end 
    
    % compute the precondition matrix 
    if nargin > 5 && strcmp(varargin{1},'none')
        C = eye(n); 
    else
        C = preconditionMatrix(df,X_,ind,n);
    end
    
    % divide variables
    [u1_,u2_] = Extract(x_,ind,n);
    [U1,~] = Extract(X,ind,n);
    [U1_,U2_] = Extract(X_,ind,n);
    [J1,J2] = Extract(C*J,ind,n);
    
    % initialization
    b = C*Y_ - C*f(x_); b = b(1:n);
    d = [inf,inf];
    
    while d(2) <= mu*d(1) && in(U1,U1_)
       
        try
            t = InvDiag(J1) * (b - OffDiag(J1)*(U1_-u1_)) - J2*(U2_ - u2_);
        catch
            break;
        end
        
        if in(U1_,u1_ + t)
            res = true;
            return; 
        end
        
        d(1) = d(2);
        d(2) = dist(U1_,u1_ + tau*t);
        U1_ = u1_ + tau*t; 
        
        temp = X_;
        temp(ind) = cartProd(U1_,U2_);
        [J1,J2] = Extract(C*df(temp),ind,n);
    end

    % try wihtout preconditioning matrix
    if any(any(C - eye(n)))
        res = containsInterval(f,df,X,X_,Y_,'none');
    else
        res = false;
    end
end

function D = InvDiag(A)
% compute the inverse of diag(A)
    D = interval(zeros(size(A)));
    for i = 1:length(A)
       D(i,i) = 1./A(i,i); 
    end
end

function A = OffDiag(A)
% compute the matrix of off-diagonal entries A - diag(A)
    for i = 1:length(A)
       A(i,i) = interval(0);  
    end
end

function d = dist(int1,int2)
% compute the distance of two intervals
    a1 = abs(supremum(int1)-supremum(int2));
    a2 = abs(infimum(int1)-infimum(int2));
    temp = max([a1,a2],[],2);
    d = sqrt(sum(temp.^2));
end

function ind = getSuitableSubmatrix(J,n)
% extract a submatrix of full rank from the Jacobian matrix

    % determine intervals that do not contain 0
    nonZero = zeros(size(J));
    
    for i = 1:size(J,1)
        for j = 1:size(J,2)
            nonZero(i,j) = ~in(J(i,j),0);
        end
    end
    
    % sort so that columns with least non-zero containments are first
    [~,index] = sort(sum(nonZero,2));

    % determine a suitable reordering of the Jacobian columns so that the
    % Jacobian matrix has full rank
    ind = zeros(1,n);
    
    for i = 1:n
        for j = 1:size(J,2)
            if ~ismember(j,ind) && ~in(J(index(i),j),0)
               ind(index(i)) = j; 
            end
        end
    end
    
    ind_ = setdiff(1:size(J,2),ind);
    
    ind = [ind,ind_];
end

function C = preconditionMatrix(df,X,ind,n)
% compute preconditioning matrix according to Sec. 4.2 in [2] so that the
% modified inclusion criterion in Eq. (10) in [2] can be applied

    % get Jacobian matrix
    x_ = center(X);
    J = df(x_);
    
    % construct extended Jacobian matrix
    J_ = [J(:,ind); zeros(size(J,2)-n,n),eye(size(J,2)-n)];
    
    % compute the inverse
    if rank(J_) < size(J_,1)
        C_ = eye(size(J_,1));
    else
        C_ = inv(J_);
    end
    
    % extract submatrix
    C = C_(1:n,1:n);
end

function [x1,x2] = Extract(x,ind,n)
% divide the variables according to the ordering defined by variable "ind"
    
    if size(x,2) > 1    % matrix
        x = x(:,ind);
        x1 = x(:,1:n); x2 = x(:,n+1:end);
    else
        x = x(ind);
        x1 = x(1:n); x2 = x(n+1:end);
    end
end

function x_ = getFactorDomain(f,Y,X)
% try to find x_ \in X such that f(x_) = center(Y) using nonlinear
% programming

    options = optimoptions('fmincon','Display','off');

    objFun = @(x) sum(x.^2);
    conFun = @(x) deal([],f(x)-center(Y));
        
    x_ = fmincon(objFun,center(X),[],[],[],[],infimum(X), ...
                           supremum(X),conFun,options);
end

function list = splitLongestGen_(pZ)
% split the longest generator, where the independent generators are also
% considered

    % determine generator with maximum length
    G = [pZ.G,pZ.Grest];
    len = sum(G.^2,1);
    [~,ind] = max(len);
    
    if ind <= size(pZ.G,2)                  % split dependent generator
        
        list = splitLongestGen(pZ);
        
    else                                    % split independent generator
        
        % get object properties
        c = pZ.c;
        Grest = pZ.Grest;
        
        % split the longest independent generator
        ind = ind - size(pZ.G,2);
        
        c1 = c + 0.5*Grest(:,ind);
        c2 = c - 0.5*Grest(:,ind);
        Grest(:,ind) = 0.5*Grest(:,ind);
        
        % construct the splitted polynomial zonotopes
        list = cell(2,1);
        
        list{1} = polyZonotope(c1,pZ.G,Grest,pZ.expMat);
        list{2} = polyZonotope(c2,pZ.G,Grest,pZ.expMat);
    end
end

%------------- END OF CODE --------------