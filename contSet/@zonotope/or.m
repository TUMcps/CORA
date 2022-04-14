function Z = or(Z1, varargin)
% or - computes an over-approximation for the union of zonotopes
%
% Syntax:  
%    Z = or(Z1, Z2)
%    Z = or(Z1, ... , Zm)
%    Z = or(Z1, ... , Zm, alg)
%    Z = or(Z1, ... , Zm, alg, order)
%
% Inputs:
%    Z1,...,Zm - zonotope objects
%    alg - algorithm used to compute the union
%               - 'linProg' (default)
%               - 'tedrake'
%               - 'iterative'
%               - 'althoff'
%               - 'parallelotope'
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    Z - resulting zonotope object enclosing the union
%
% Example: 
%    zono1 = zonotope([4 2 2;1 2 0]);
%    zono2 = zonotope([3 1 -1 1;3 1 2 0]);
%
%    res = zono1 | zono2;
%
%    figure
%    hold on
%    plot(zono1,[1,2],'r');
%    plot(zono2,[1,2],'b');
%    plot(res,[1,2],'g');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      20-September-2013
% Last update:  12-November-2019 (NK: added two new algorithms)
% Last revision:---

%------------- BEGIN CODE --------------

    % default values
    alg = 'linProg';
    order = [];

    % distinguish case with two and case with more sets
    if nargin == 2 || (nargin > 2 && ...
                      (isempty(varargin{2}) || ischar(varargin{2})))
                  
        Z2 = varargin{1};

        % determine zonotope object
        if ~isa(Z1,'zonotope')
            temp = Z1;
            Z1 = Z2;
            Z2 = temp;
        end

        % different cases depending on the class of the second set
        if isa(Z2,'zonotope') || isa(Z2,'interval') || isnumeric(Z2)
            
            if ~isa(Z2,'zonotope')
               Z2 = zonotope(Z2); 
            end
            
            % parse input arguments
            if nargin > 2 && ~isempty(varargin{2})
                alg = varargin{2};
            end
            
            if nargin > 3 && ~isempty(varargin{3})
                order = varargin{3};
            end
            
            % compute over-approximation of the union with the selected 
            % algorithm
            if strcmp(alg,'linProg')
               Z = unionLinProg({Z1,Z2},order);
            elseif strcmp(alg,'althoff')
               Z = unionAlthoff(Z1,{Z2},order);
            elseif strcmp(alg,'tedrake')
               Z = unionTedrake({Z1,Z2},order);
            elseif strcmp(alg,'iterative')
               Z = unionIterative({Z1,Z2},order);
            elseif strcmp(alg,'parallelotope')
               Z = unionParallelotope({Z1,Z2});
            else
               error('Wrong value for input argument ''alg''!'); 
            end

        elseif isa(Z2,'mptPolytope') || isa(Z2,'conZonotope') || ...
               isa(Z2,'zonoBundle') || isa(Z2,'conPolyZono')

            Z = Z2 | Z1;     

        else
            % throw error for given arguments
            error(noops(Z1,Z2));
        end
    
    else
        
        % parse input arguments
        Zcell = {};
        counter = [];

        for i = 2:nargin
           if isa(varargin{i-1},'zonotope')
              Zcell{end+1,1} = varargin{i-1}; 
           else
              counter = i;
              break;
           end
        end

        if ~isempty(counter)
            if nargin >= counter && ~isempty(varargin{counter-1})
                alg = varargin{counter-1}; 
            end
            if nargin >= counter+1 && ~isempty(varargin{counter})
                order = varargin{counter};
            end
        end

        % compute over-approximation of the union with the selected algorithm
        if strcmp(alg,'linProg')
           Z = unionLinProg([{Z1};Zcell],order);
        elseif strcmp(alg,'althoff')
           Z = unionAlthoff(Z1,Zcell,order);
        elseif strcmp(alg,'tedrake')
           Z = unionTedrake([{Z1};Zcell],order);
        elseif strcmp(alg,'iterative')
           Z = unionIterative([{Z1};Zcell],order);
        elseif strcmp(alg,'parallelotope')
           Z = unionParallelotope([{Z1};Zcell]);
        else
           error('Wrong value for input argument ''alg''!'); 
        end
    end
end



% Auxiliary Functions -----------------------------------------------------

function Z = unionTedrake(Zcell,order)
% compute the union by solving a linear program with respect to the
% zonotope containment constraints presented in [1]

    % construct generator matrix of the final zonotope
    Z_ = unionIterative(Zcell,order);
    G = generators(Z_);
    
    Y = G*diag(1./sqrt(sum(G.^2,1)));
    [n,ny] = size(Y);
    Hy = [eye(ny);-eye(ny)];

    % construct linear constraints for each zonotope
    Aeq = [];
    beq = [];
    A = [];
    A_ = [];
    
    for i = 1:length(Zcell)
       
        % obtain generator matrix and center from the current zonotope
        X = generators(Zcell{i});
        x = center(Zcell{i});
        nx = size(X,2);
        Hx = [eye(nx);-eye(nx)];
        hx = ones(2*nx,1);
        
        % construct constraint X = Y * T
        temp = repmat({Y},[1,nx]);
        A1 = [blkdiag(temp{:}),zeros(n*nx,4*nx*ny),zeros(n*nx,ny)];
        b1 = reshape(X,[n*nx,1]);

        % construct constraint x = Y * beta
        A2 = [zeros(n,nx*ny),zeros(n,4*nx*ny),Y];
        b2 = -x;

        % construct constraint lambda * Hx = Hy * T
        Atemp = [];
        for j = 1:size(Hx,2)
            h = Hx(:,j);
            temp = repmat({h'},[1,size(Hy,1)]);
            Atemp = [Atemp;blkdiag(temp{:})];
        end
        
        temp = repmat({Hy},[1,size(Hx,2)]);
        A3 = [blkdiag(temp{:}),-Atemp,zeros(size(Atemp,1),ny)];
        b3 = zeros(size(A3,1),1);
        
        % add current equality constraint to overall equality constraints
        Atemp = [A1;A2;A3];
        btemp = [b1;b2;b3];
        
        Aeq = blkdiag(Aeq,Atemp);
        beq = [beq;btemp];
        
        % construct constraint lambda * hx <= hy + Hy beta
        temp = repmat({hx'},[1,size(Hy,1)]);
        A_ = [A_;-eye(size(Hy,1))];
        A1 = [zeros(size(Hy,1),size(Y,2)*size(X,2)),blkdiag(temp{:}),-Hy];
        
        % construct constraint lambda >= 0
        A2 = [zeros(4*nx*ny,nx*ny),-eye(4*nx*ny),zeros(4*nx*ny,ny)];
        A_ = [A_;zeros(4*nx*ny,2*ny)];
        
        A = blkdiag(A,[A1;A2]);
    end
    
    % solve linear program
    f = [ones(2*ny,1);zeros(size(Aeq,2),1)];
    
    Atemp = [-eye(ny),-eye(ny)];
    A = [[A_,A];[Atemp,zeros(ny,size(A,2))]];
    b = zeros(size(A,1),1);
    Aeq = [zeros(size(Aeq,1),2*ny),Aeq];
    
    val = linprog(f',A,b,Aeq,beq,[],[]);
    
    % construct the resulting zonotope
    ub = val(1:ny);
    lb = -val(ny+1:2*ny);
    int = interval(lb,ub);
    
    c = Y*center(int);
    G = Y * diag(rad(int));
    
    Z = zonotope([c,G]);

end

function Z = unionLinProg(Zcell,order)
% compute an enclosing zonotope using linear programming. As the
% constraints for the linear program we compute the upper and lower bound
% for all zonotopes that are enclosed in the normal directions of the
% halfspace representation of the enclosing zonotope

    % construct generator matrix of the final zonotope
    Z_ = unionIterative(Zcell,order);
    G = generators(Z_);
    
    G = G*diag(1./sqrt(sum(G.^2,1)));
    
    % compute the directions of the boundary halfspaces
    [n,m] = size(G);
    
    Z = zonotope([zeros(n,1),G]);
    Z = halfspace(Z);
    
    C = Z.halfspace.H;
    
    % compute bounds for each halfspace
    d = zeros(size(C,1),1);
    
    for i = 1:size(C,1)
       
        val = -inf;
        
        % loop over all zonotopes
        for j = 1:length(Zcell)
           
            % compute bound for the current zonotope
            Ztemp = C(i,:)*Zcell{j};
            valTemp = Ztemp.Z(1) + sum(abs(Ztemp.Z(2:end)));
            
            % update bound
            val = max(val,valTemp);
        end
        
        d(i) = val;
    end

    % solve linear program
    f = [zeros(n,1);ones(m,1)];
    
    A = [];
    b = -d;
    
    for i = 1:size(C,1)
       A = [A;-[C(i,:),abs(C(i,:)*G)]];
    end
    
    A = [A;[zeros(m,n),-eye(m)]];
    b = [b;zeros(m,1)];
    
    options = optimoptions('linprog','display','off');
    
    x = linprog(f',A,b,[],[],[],[],options);
    
    % construct final zonotope
    c = x(1:n);
    scal = x(n+1:end);
    
    Z = zonotope([c,G*diag(scal)]);
      
end

function Z = unionIterative(Zcell,order)
% iteratively unite all zonotopes with the enclose function
    
    Zcell_ = cell(length(Zcell),1);
    
    % loop until all sets have been united
    while length(Zcell) > 1
       
        counter = 1;
        counter_ = 1;
        
        % loop over all sets in the current queue
        while counter < length(Zcell)
           
            % unite sets using function enclose
            Zcell_{counter_} = enclose(Zcell{counter},Zcell{counter+1});
            
            counter_ = counter_ + 1;
            counter = counter + 2;           
        end
        
        % consider remaining set
        if counter == length(Zcell)
           Zcell_{counter_} = Zcell{counter}; 
           counter_ = counter_ + 1;
        end
        
        Zcell = Zcell_(1:counter_-1);
    end
    
    % reduce the resulting zonotope to the desired order
    Z = Zcell{1};
    
    if ~isempty(order)
        if order == 1
            Z = reduce(Z,'methC',order);
        else
            Z = reduce(Z,'pca',order);
        end
    end
    
end

function Z = unionAlthoff(Z1,Zcell,order)

    % init
    Zmat = [];

    % dimension
    n = length(center(Z1));

    % obtain minimum number of generators every zonotope has
    minNrOfVecs = length(Z1.Z(1,:));
    for iSet = 1:length(Zcell)
        minNrOfVecs = min(minNrOfVecs, length(Zcell{iSet}.Z(1,:)));
    end

    % obtain Zcut
    Zcut{1} = Z1.Z(:, 1 : minNrOfVecs);
    for iSet = 1:length(Zcell)
        Zcut{iSet+1} = Zcell{iSet}.Z(:, 1 : minNrOfVecs);
    end

    % obtain Zadd
    Zadd{1} = Z1.Z(:, minNrOfVecs+1 : end);
    for iSet = 1:length(Zcell)
        Zadd{iSet+1} = Zcell{iSet}.Z(:, minNrOfVecs+1 : end);
    end

    % compute vertex sets for each set of generators
    for iGen = 1:minNrOfVecs
        v(:,1) = Zcut{1}(:,iGen);
        for iSet = 2:length(Zcut)
            v_new = Zcut{iSet}(:,iGen);
            % if direction is correct
            if v(:,1).'*v_new > 0
                v(:, iSet) = v_new;
            % flip direction
            else
                v(:, iSet) = - v_new;
            end
        end

        % compute vertices
        V = vertices(v);

        % compute enclosing zonotope
        Z_encl = zonotope.enclosePoints(V);

        % concatenate enclosing zonotopes
        Zmat(:,end+1 : end+n+1) = Z_encl.Z;
    end

    % add Zadd to the resulting generator matrix
    for i = 1:length(Zadd)
       Zmat = [Zmat,Zadd{i}]; 
    end

    % create enclosing zonotope
    Z = zonotope(Zmat);
    
    % reduce zonotope to the desired zonotope order
    if ~isempty(order)
       Z = reduce(Z,'pca',order); 
    end
end

function Z = unionParallelotope(Zcell)
% enclose the union of the zonotopes with a parallelotope

    % obtain matrix of points from generator matrices
    V = [];
    
    for i = 1:length(Zcell)
        G = generators(Zcell{i});
        V = [V,G,-G];
    end

    % compute the arithmetic mean of the vertices
    mean = sum(V,2)/length(V(1,:));

    % obtain sampling matrix
    translation = mean*ones(1,length(V(1,:)));
    sampleMatrix = V-translation;

    % compute the covariance matrix
    C = cov(sampleMatrix');

    % singular value decomposition
    [U,~,~] = svd(C);

    % enclose zonotopes with intervals in the transformed space
    int = [];
    
    for i = 1:length(Zcell)
        temp = U.'*Zcell{i};
        int = interval(temp) | int;
    end

    % transform back to original space
    Z = U*zonotope(int);
end

%------------- END OF CODE --------------