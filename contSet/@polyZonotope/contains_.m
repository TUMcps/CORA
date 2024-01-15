function res = contains_(pZ,S,type,varargin)
% contains_ - determines if a polynomial zonotope contains a set or a point
%
% Syntax:
%    res = contains_(pZ,S)
%    res = contains_(pZ,S,type)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object or single point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);
%
%    p1 = [1;1];
%    p2 = [-1;3];
%    pZ1 = polyZonotope([0;0],0.3*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],0.4*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
% 
%    contains(pZ,p1,'approx')
%    contains(pZ,p2,'approx')
%    contains(pZ,pZ1,'approx')
%    contains(pZ,pZ2,'approx')
% 
%    figure; hold on;
%    plot(pZ,[1,2],'b');
%    plot(p1(1),p1(2),'.g','MarkerSize',20);
%    plot(p2(1),p2(2),'.r','MarkerSize',20);
%    plot(pZ1,[1,2],'g');
%    plot(pZ2,[1,2],'r');
%
% References: 
%   [1] O. Mullier and et al. "General Inner Approximation of Vector-valued 
%       Functions", Reliable Computing, 2013
%   [2] N. Kochdumper "Extensions of Polynomial Zonotopesand their 
%       Application to Verification of Cyber-Physical Systems", 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains, conZonotope/contains

% Authors:       Niklas Kochdumper
% Written:       13-January-2020 
% Last update:   25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

    % init result
    res = false;

    % check user inputs 
    if strcmp(type,'exact')
        throw(CORAerror('CORA:noExactAlg',pZ,S));
    end
    
    % initialize variables
    fHan = @(x) aux_funcPoly(x,pZ.c,pZ.G,pZ.GI,pZ.E);
    jacHan = aux_funHanJacobian(pZ.G,pZ.GI,pZ.E);
        
    temp = ones(length(pZ.id) + size(pZ.GI,2),1);
    X = interval(-temp,temp);
        
    % point in polynomial zonotope containment
    if isnumeric(S)
        
        % init containment check
        res = false(1,size(S,2));
        for i=1:size(S,2)
            % quick check: point is center of polyZonotope?
            res(i) = all(withinTol(S(:,i),pZ.c));

            if ~res(i)
                % try to prove that the point is inside the polynomial zonotope
                % using the approach from [1]
                Y = interval(S(:,i));
                
                x = aux_getFactorDomain(fHan,Y,X);
                X_ = interval(x - 1e-10*ones(size(x)), x + 1e-10*ones(size(x)));
                
                res(i) = aux_containsInterval(fHan,jacHan,X,X_,Y);
            end
        end
        
    % set in zonotope containment
    elseif isa(S,'interval') || isa(S,'zonotope') || ...
           isa(S,'polytope') || isa(S,'zonoBundle') || ...
           isa(S,'conZonotope') || isa(S,'taylm') || ...
           isa(S,'polyZonotope') || isa(S,'conPolyZono')
        
        % convert set to polynomial zonotope
        if ~isa(S,'polyZonotope')
            S = polyZonotope(S);
        end
        
        % try to prove that set containment does not hold using 
        % Proposition 3.1.34 in [2]
        if aux_disproveContainment(pZ,S)
            return; 
        end
        
        % try to prove that set containment holds using 
        % Proposition 3.1.36 in [2]
        res = aux_proveContainment(S,fHan,jacHan,X);

    else
        throw(CORAerror('CORA:noops',pZ,S));
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_disproveContainment(pZ1,pZ2)
% use contraction to prove that pZ2 is not a subset of pZ1 
% (see Proposition 3.1.34 in [2])

    % construct polynomial constraint for the intersection
    c = pZ1.c - pZ2.c;
    G = [pZ1.G -pZ2.G];
    GI = [pZ1.GI,-pZ2.GI];
    E = blkdiag(pZ1.E,pZ2.E);
        
    % contract the factor domain \alpha_k in [-1,1] based on the
    % polynomial constraint
    n = size(E,1) + size(GI,2);
    dom = interval(-ones(n,1),ones(n,1));

    I = contractPoly(c,G,GI,E,dom,'all',3);

    % pZ2 is proven to be not contained in pZ1 -> res = true
    ind = size(pZ1.E,1) + [(1:size(pZ2.E,1))'; ...
                     size(pZ2.E,1) + size(pZ1.GI,2) + (1:size(pZ2.GI,2))']; 
    temp = I(ind);
    res = false;
    
    if representsa_(temp,'emptySet',eps) || any(supremum(temp) < 1) || any(infimum(temp) > -1)
       res = true; 
    end
end

function res = aux_proveContainment(obj,fHan,jacHan,X)
% try to prove that the set is contained in the polynomial zonotope by
% recursiveley splitting the set and using the method in [1] to prove that
% an interval enclosure of the splitted sets is contained in the polynomial
% zonotope (see Proposition 3.1.36 in [2])

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
        
            x = aux_getFactorDomain(fHan,Y,X);
            X_ = x + pinv(jacHan(x))*(Y - center(Y));
        
            temp = aux_containsInterval(fHan,jacHan,X,X_,Y); 
            
            % split the set if it is not contained
            if ~temp
               sets = aux_splitLongestGen_(list{j});
               
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

function f = aux_funcPoly(x,c,G,GI,E)
% evaluate the mapping function of a polynomial zonotope

    % initialization
    f = c;
    n = size(E,1);
    x1 = x(1:n); x2 = x(n+1:end);
    
    % loop over all dependent generators
    for i = 1:size(G,2)
        f = f + G(:,i)*prod(x1.^E(:,i),1); 
    end

    % add independent generators
    if ~isempty(GI)
        f = f + GI*x2;
    end
end

function J = aux_jacobianPoly(x,G,GI,E,n,p)
% compute the jacobian matrix for the mapping function of a polynomial
% zonotope

    % initialization
    J = 0 * repmat(x(1),[n,length(E)]);
    x_ = x(1:p);
    
    % loop over all variables
    for i = 1:length(E)
       for j = 1:size(E{i},2)
          J(:,i) = J(:,i) + G{i}(:,j)*prod(x_.^E{i}(:,j),1); 
       end
    end
    
    % consider independent generators
    J = [J,GI];
end

function jacHan = aux_funHanJacobian(G,GI,E)
% compute a function handle for the jacobian matrix for the mapping 
% function of a polynomial zonotope

    % compute exponent matrix differentiazed for each variable
    Elist = cell(size(E,1),1);
    Glist = cell(size(E,1),1);

    for i = 1:length(Elist)
       ind = find(E(i,:) > 0);
       temp = E(:,ind);
       temp(i,:) = temp(i,:) - 1;
       Elist{i} = temp;
       Glist{i} = G(:,ind) * diag(E(i,ind));
    end

    % function handle for jacobian matrix
    jacHan=@(x)aux_jacobianPoly(x,Glist,GI,Elist,size(G,1),size(E,1));
end

function res = aux_containsInterval(f,df,X,X_,Y_,varargin)
% implementation of the function Inner() according to Alg. 3 in [2] which
% checks if an interval is part of the image of a nonlinear function

    % parameters
    tau = 1.01;
    mu = 0.9;
    
    % dermine suitable dimensions of the domain X
    J = df(X_);
    x_ = center(X_);
    n = size(J,1);
    
    ind = aux_getSuitableSubmatrix(J,n);
    
    % check if Jacobian matrix has full rank
    if ismember(0,ind)  
        res = 0;
        return;
    end 
    
    % compute the precondition matrix 
    if nargin > 5 && strcmp(varargin{1},'none')
        C = eye(n); 
    else
        C = aux_preconditionMatrix(df,X_,ind,n);
    end
    
    % divide variables
    [u1_,u2_] = aux_Extract(x_,ind,n);
    [U1,~] = aux_Extract(X,ind,n);
    [U1_,U2_] = aux_Extract(X_,ind,n);
    [J1,J2] = aux_Extract(C*J,ind,n);
    
    % initialization
    b = C*Y_ - C*f(x_); b = b(1:n);
    d = [inf,inf];
    
    while d(2) <= mu*d(1) && contains_(U1,U1_,'exact',0)
       
        try
            t = aux_InvDiag(J1) * (b - aux_OffDiag(J1)*(U1_-u1_)) - J2*(U2_ - u2_);
        catch
            break;
        end
        
        if contains_(U1_,u1_ + t,'exact',0)
            res = true;
            return; 
        end
        
        d(1) = d(2);
        d(2) = aux_dist(U1_,u1_ + tau*t);
        U1_ = u1_ + tau*t; 
        
        temp = X_;
        temp(ind) = cartProd_(U1_,U2_,'exact');
        [J1,J2] = aux_Extract(C*df(temp),ind,n);
    end

    % try without pre-conditioning matrix
    if any(any(C - eye(n)))
        res = aux_containsInterval(f,df,X,X_,Y_,'none');
    else
        res = false;
    end
end

function D = aux_InvDiag(A)
% compute the inverse of diag(A)
    D = interval(zeros(size(A)));
    for i = 1:length(A)
       D(i,i) = 1./A(i,i); 
    end
end

function A = aux_OffDiag(A)
% compute the matrix of off-diagonal entries A - diag(A)
    for i = 1:length(A)
       A(i,i) = interval(0);  
    end
end

function d = aux_dist(int1,int2)
% compute the distance of two intervals
    a1 = abs(supremum(int1)-supremum(int2));
    a2 = abs(infimum(int1)-infimum(int2));
    temp = max([a1,a2],[],2);
    d = sqrt(sum(temp.^2));
end

function ind = aux_getSuitableSubmatrix(J,n)
% extract a submatrix of full rank from the Jacobian matrix

    % determine intervals that do not contain 0
    nonZero = zeros(size(J));
    
    for i = 1:size(J,1)
        for j = 1:size(J,2)
            nonZero(i,j) = ~contains(J(i,j),0);
        end
    end
    
    % sort so that columns with least non-zero containments are first
    [~,index] = sort(sum(nonZero,2));

    % determine a suitable reordering of the Jacobian columns so that the
    % Jacobian matrix has full rank
    ind = zeros(1,n);
    
    for i = 1:n
        for j = 1:size(J,2)
            if ~ismember(j,ind) && ~contains(J(index(i),j),0)
               ind(index(i)) = j; 
            end
        end
    end
    
    ind_ = setdiff(1:size(J,2),ind);
    
    ind = [ind,ind_];
end

function C = aux_preconditionMatrix(df,X,ind,n)
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

function [x1,x2] = aux_Extract(x,ind,n)
% divide the variables according to the ordering defined by variable "ind"
    
    if size(x,2) > 1    % matrix
        x = x(:,ind);
        x1 = x(:,1:n); x2 = x(:,n+1:end);
    else
        x = x(ind);
        x1 = x(1:n); x2 = x(n+1:end);
    end
end

function x_ = aux_getFactorDomain(f,Y,X)
% try to find x_ \in X such that f(x_) = center(Y) using nonlinear
% programming

    options = optimoptions('fmincon','Display','off');

    objFun = @(x) sum(x.^2);
    conFun = @(x) deal([],f(x)-center(Y));
        
    x_ = fmincon(objFun,center(X),[],[],[],[],infimum(X), ...
                           supremum(X),conFun,options);
end

function list = aux_splitLongestGen_(pZ)
% split the longest generator, where the independent generators are also
% considered

    % determine generator with maximum length
    G = [pZ.G,pZ.GI];
    len = sum(G.^2,1);
    [~,ind] = max(len);
    
    if ind <= size(pZ.G,2)                  % split dependent generator
        
        list = splitLongestGen(pZ);
        
    else                                    % split independent generator
        
        % get object properties
        c = pZ.c;
        GI = pZ.GI;
        
        % split the longest independent generator
        ind = ind - size(pZ.G,2);
        
        c1 = c + 0.5*GI(:,ind);
        c2 = c - 0.5*GI(:,ind);
        GI(:,ind) = 0.5*GI(:,ind);
        
        % construct the splitted polynomial zonotopes
        list = cell(2,1);
        
        list{1} = polyZonotope(c1,pZ.G,GI,pZ.E);
        list{2} = polyZonotope(c2,pZ.G,GI,pZ.E);
    end
end

% ------------------------------ END OF CODE ------------------------------
