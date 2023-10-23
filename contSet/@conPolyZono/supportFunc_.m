function val = supportFunc_(cPZ,dir,type,method,splits,varargin)
% supportFunc_ - calculate the upper or lower bound of a constrained 
%    polynomial zonotope object along a certain direction
%
% Syntax:
%    val = supportFunc_(cPZ,dir)
%    val = supportFunc_(cPZ,dir,type)
%    val = supportFunc_(cPZ,dir,type,method)
%    val = supportFunc_(cPZ,dir,type,'split',splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    method - method that is used to calculate the bounds for the dependent
%             part of the constrained polynomial zonotope
%              'conZonotope': conversion to a constrained zonotope
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'quadProg': quadratic programming
%    split - number of splits that are performed to calculate the bounds
%
% Outputs:
%    val - interval object specifying the upper and lower bound along the
%          direction
%
% Examples:
%    c = [0;0];
%    G = [2 0; 0 1];
%    E = [1 0; 0 1];
%    A = [1 1];
%    b = 1;
%    EC = [2 0; 0 2];
%    GI = [1;0];
%    cPZ = conPolyZono(c,G,E,A,b,EC,GI);
%
%    dir = [1;1];
%    
%    b1 = supportFunc(cPZ,dir,'range');
%    ub1 = conHyperplane(dir,supremum(b1));
%    lb1 = conHyperplane(dir,infimum(b1));
%
%    b2 = supportFunc(cPZ,dir,'range','split');
%    ub2 = conHyperplane(dir,supremum(b2));
%    lb2 = conHyperplane(dir,infimum(b2));
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(ub1,[1,2],'b');
%    plot(lb1,[1,2],'b');
%    plot(ub2,[1,2],'g');
%    plot(lb2,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, polyZonotope/supportFunc_, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper
% Written:       29-July-2018
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------
    
    % different methods to calculate the over-approximation
    if strcmp(method,'interval')
       
        % project the polynomial zonotope onto the direction
        cPZ_ = dir' * cPZ;
        
        % interval enclosure
        val = interval(zonotope(cPZ_));
        
        if strcmp(type,'lower')
            val = infimum(val); 
        elseif strcmp(type,'upper')
            val = supremum(val);
        elseif strcmp(type,'range')
            % val = val
        end
        
    elseif strcmp(method,'split')
        
        val = aux_supportFuncSplit(cPZ,dir,type,splits);
        
    elseif strcmp(method,'conZonotope')
        
        val = aux_supportFuncConZonotope(cPZ,dir,type);
                
    elseif strcmp(method,'quadProg')
        
        val = aux_supportFuncQuadProg(cPZ,dir,type);
        
    end
end


% Auxiliary functions -----------------------------------------------------

function val = aux_supportFuncSplit(cPZ,dir,type,splits)
% comptue the support function by recursively splitting the set

    % project the polynomial zonotope onto the direction
    cPZ_ = dir' * cPZ;

    % split the constrained polynomial zonotope multiple times to 
    % obtain a better over-approximation of the real shape
    cPZsplit{1} = cPZ_;

    for i = 1:splits
        qZnew = [];
        for j = 1:length(cPZsplit)
            res = splitLongestGen(cPZsplit{j});
            qZnew{end+1} = res{1};
            qZnew{end+1} = res{2};
        end
        cPZsplit = qZnew;
    end

    % calculate the interval over-approximation
    Min = Inf;
    Max = -Inf;

    for i = 1:length(cPZsplit)
       try
           I = interval(cPZsplit{i},'interval');
       catch ME
           if strcmp(ME.identifier,'CORA:emptySet')
              continue;
           else
              rethrow(ME);
           end
       end
       Min = min(Min,infimum(I));
       Max = max(Max,supremum(I));
    end

    if isinf(Min) || isinf(Max)
        throw(CORAerror('CORA:emptySet'));
    end

    % extract the desired bound
    if strcmp(type,'lower')
       val = Min; 
    elseif strcmp(type,'upper')
       val = Max;
    elseif strcmp(type,'range')
       val = interval(Min,Max);
    end
end

function val = aux_supportFuncConZonotope(cPZ,dir,type)
% compute the support function based on a conZonotope enclosure

    % project the polynomial zonotope onto the direction
    cPZ_ = dir' * cPZ;

    % split into dependent part and independent part
    GI = cPZ_.GI;
    cPZ_.GI = [];

    % enlose dependent part with a constrained zonotope
    cZ = conZonotope(cPZ_);

    % compute support function of the constrained zonotope
    if strcmp(type,'range')
        l = supportFunc_(cZ,1,'lower');
        u = supportFunc_(cZ,1,'upper');
        val1 = interval(l,u);
    else
        val1 = supportFunc_(cZ,1,type);
    end

    % compute support function of the independent part
    if ~isempty(GI)
        Z = zonotope(0,GI);
        if strcmp(type,'range')
            l = supportFunc_(Z,1,'lower');
            u = supportFunc_(Z,1,'upper');
            val2 = interval(l,u);
        else
            val2 = supportFunc_(Z,1,type);
        end
    else
        val2 = 0*val1;
    end

    % combine the results
    val = val1 + val2;
end

function val = aux_supportFuncQuadProg(cPZ,dir,type)
% compute the support function based on quadratic programming

    % consider different types of bounds
    if strcmp(type,'range')
        l = supportFunc_(cPZ,dir,'lower','quadProg',[]);
        u = supportFunc_(cPZ,dir,'upper','quadProg',[]);
        val = interval(l,u);
        return;
    end

    if isempty(cPZ.A)
       pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E);
       val = supportFunc_(pZ,dir,type,'quadProg',[],[]);
       return;
    end
    
    % project the polynomial zonotope onto the direction
    cPZ_ = dir' * cPZ;
    
    c = cPZ_.c; cPZ_.c = 0;
    if strcmp(type,'upper')
        cPZ_ = (-1)*cPZ_;
    end

    % extract linear part f*x for states
    indLin = find(sum(cPZ_.E,1) <= 1);
    p = length(cPZ_.id); f = zeros(p,1);

    for i = 1:length(indLin)
        ind = find(cPZ_.E(:,indLin(i)) == 1);
        f(ind(1)) = cPZ_.G(:,indLin(i));
    end

    % extract quadratic part x'*H*x for states
    indQuad = find(sum(cPZ_.E,1) == 2);
    H = zeros(p);

    for i = 1:length(indQuad)
       if max(cPZ_.E(:,indQuad(i))) == 2
           ind = find(cPZ_.E(:,indQuad(i)) == 2);
           H(ind(1),ind(1)) = cPZ_.G(:,indQuad(i));
       else
           ind = find(cPZ_.E(:,indQuad(i)) == 1);
           H(ind(1),ind(2)) = 0.5 * cPZ_.G(:,indQuad(i));
           H(ind(2),ind(1)) = 0.5 * cPZ_.G(:,indQuad(i));
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
        val = supportFunc_(cPZ,dir,type,'conZonotope',[],[]);
        return;
    end
    if ~isempty(ind_)
        temp = zeros(size(d)); temp(ind_) = d(ind_);
        N = N + V*diag(temp)*V';
    end

    % extract linear part A*x = b for constraints
    indLinCon = find(sum(cPZ_.EC,1) <= 1);
    A = zeros(size(cPZ_.A,1),p);

    for i = 1:length(indLinCon)
        ind = find(cPZ_.EC(:,indLinCon(i)) == 1);
        A(:,ind(1)) = cPZ_.A(:,indLinCon(i));
    end

    % enclose remaining part with additional factors
    ind = setdiff(1:size(cPZ_.G,2),[indQuad,indLin]);
    indCon = setdiff(1:size(cPZ_.A,2),indLinCon);

    G = cPZ_.G(:,ind); E = cPZ_.E(:,ind);
    A_ = cPZ_.A(:,indCon); EC = cPZ.EC(:,indCon);

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
    if isempty(E)
       if ~isempty(EC)
           Z = zonotope(polyZonotope(-cPZ_.b,A_,[],EC));

           H = blkdiag(H,zeros(size(Z.G,2)));
           f = [f; zeros(size(Z.G,2),1)];
           A = [A, Z.G]; b = -Z.c;
       else
           b = cPZ_.b;
       end
    else
        if ~isempty(EC)
            cZ = conZonotope(conPolyZono(c,G,E,A_,cPZ_.b, ...
                             EC),'extend','none');

            H = blkdiag(H,zeros(size(cZ.G,2)));
            f = [f; cZ.G']; c = cZ.c;
            A = [A, cZ.A]; b = cZ.b;
        else
            Z = zonotope(polyZonotope(c,G,[],E));

            H = blkdiag(H,zeros(size(Z.G,2)));
            f = [f; Z.G']; c = cPZ_.c;
            A = [A, zeros(size(A,1),size(Z.G,2))]; b = cPZ.b;
        end
    end

    % solve quadratic program
    options = optimoptions('quadprog','Display','off');
    w = warning(); warning('off');

    [~,val] = quadprog(2*H,f,[],[],A,b,-ones(length(f),1), ...
                       ones(length(f),1),[],options);

    warning(w);

    val = val + c;

    if strcmp(type,'upper')
       val = -val; 
    end

    % compute bound for independent generators
    if ~isempty(cPZ.GI)
       val = val + supportFunc_( ...
           zonotope(zeros(dim(cPZ),1),cPZ.GI), dir,type); 
    end
end

% ------------------------------ END OF CODE ------------------------------
