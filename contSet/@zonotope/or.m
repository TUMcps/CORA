function S_out = or(Z,S,varargin)
% or - computes an over-approximation for the union of zonotopes
%
% Syntax:
%    S_out = Z | S;
%    S_out = or(Z,S)
%    S_out = or(Z,S,alg,order)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object, numeric, cell-array
%    alg - algorithm used to compute the union
%               - 'linprog' (default)
%               - 'tedrake'
%               - 'iterative'
%               - 'althoff'
%               - 'parallelotope'
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    S_out - zonotope object enclosing the union
%
% Example: 
%    Z1 = zonotope([4 2 2;1 2 0]);
%    Z2 = zonotope([3 1 -1 1;3 1 2 0]);
%    S_out = Z1 | Z2;
%
%    figure; hold on;
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'b');
%    plot(S_out,[1,2],'g');
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

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       20-September-2013
% Last update:   12-November-2019 (NK, added two new algorithms)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TODO: solve cell-array case...

% default values
[alg,order] = setDefaultValues({'linprog',[]},varargin);

% check input arguments
inputArgsCheck({{Z,'att',{'zonotope','numeric'}},...
                {S,'att',{'contSet','numeric'}},...
                {alg,'str',{'linprog','tedrake','iterative','althoff','parallelotope'}},...
                {order,'att','numeric'}});

% ensure that numeric is second input argument
[Z,S] = reorderNumeric(Z,S);

% check dimensions
equalDimCheck(Z,S);

% write all into a cell-array
if iscell(S)
    S_cell = [{Z}; S];
else
    S_cell = {Z; S};
end

% only zonotopes, intervals, and numeric allowed
if ~all(cellfun(@(S) isa(S,'zonotope') || isa(S,'interval') || isnumeric(S),...
        S_cell,'UniformOutput',true))
    throw(CORAerror('CORA:noops',Z,S));
end

% convert all sets to zonotopes, skip empty sets
S_cell = cellfun(@(S) zonotope(S),S_cell,'UniformOutput',false);
S_cell = S_cell(cellfun(@(S) ~representsa_(S,'emptySet',eps),S_cell,'UniformOutput',true));

% if only one set remains, this is the union
if isempty(S_cell)
    S_out = emptySet(dim(Z));
    return
elseif numel(S_cell) == 1
    S_out = S_cell{1};
    return
end

% compute over-approximation of the union with the selected algorithm
switch alg
    case 'linprog'
        S_out = aux_unionLinprog(S_cell,order);
    case 'althoff'
        S_out = aux_unionAlthoff(S_cell{1},S_cell(2:end),order);
    case 'tedrake'
        S_out = aux_unionTedrake(S_cell,order);
    case 'iterative'
        S_out = aux_unionIterative(S_cell,order);
    case 'parallelotope'
        S_out = aux_unionParallelotope(S_cell);
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_unionTedrake(Zcell,order)
% compute the union by solving a linear program with respect to the
% zonotope containment constraints presented in [1]

    % construct generator matrix of the final zonotope
    Z_ = aux_unionIterative(Zcell,order);
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

    problem.f = f';
    problem.Aineq = A;
    problem.bineq = b;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = [];
    problem.ub = [];
    
    val = CORAlinprog(problem);
    
    % construct the resulting zonotope
    ub = val(1:ny);
    lb = -val(ny+1:2*ny);
    int = interval(lb,ub);
    
    c = Y*center(int);
    G = Y * diag(rad(int));
    
    Z = zonotope([c,G]);

end

function Z = aux_unionLinprog(Zcell,order)
% compute an enclosing zonotope using linear programming. As the
% constraints for the linear program we compute the upper and lower bound
% for all zonotopes that are enclosed in the normal directions of the
% halfspace representation of the enclosing zonotope

    % construct generator matrix of the final zonotope
    nrZ = length(Zcell);
    Z_ = aux_unionIterative(Zcell,order);
    Z_ = compact_(Z_,'zeros',eps);
    G = generators(Z_);
    
    G = G*diag(1./sqrt(sum(G.^2,1)));
    
    % compute the directions of the boundary halfspaces
    [n,nrGen] = size(G);
    P = polytope(Z_ - Z_.c);
    nrIneq = length(P.b);
    
    val = zeros(nrIneq,nrZ);
    for i = 1:nrIneq
        % loop over all zonotopes
        for j = 1:nrZ
            % compute bound for the current zonotope (note: this is a
            % direct implementation of zonotope/supportFunc_ ...)
            Z_proj = P.A(i,:)*Zcell{j};
            val(i,j) = Z_proj.c + sum(abs(Z_proj.G));
        end
    end
    d = max(val,[],2);

    % solve linear program
    f = [zeros(n,1);ones(nrGen,1)];
    
    problem.f = f';
    problem.Aineq = [[-P.A, -abs(P.A*G)];[zeros(nrGen,n),-eye(nrGen)]];
    problem.bineq = [-d;zeros(nrGen,1)];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = [];
    problem.ub = [];
    
    x = CORAlinprog(problem);
    
    % construct final zonotope
    c = x(1:n);
    scal = x(n+1:end);
    
    Z = zonotope(c,G*diag(scal));
      
end

function Z = aux_unionIterative(Zcell,order)
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

function Z = aux_unionAlthoff(Z1,Zcell,order)

    % init
    Zmat = [];

    % dimension
    n = dim(Z1);

    % obtain minimum number of generators every zonotope has
    minNrOfGens = size(Z1.G,2);
    for iSet = 1:length(Zcell)
        minNrOfGens = min(minNrOfGens, size(Zcell{iSet}.G,2));
    end

    % obtain Zcut
    Zcut{1} = Z1.G(:, 1:minNrOfGens);
    for iSet = 1:length(Zcell)
        Zcut{iSet+1} = [Zcell{iSet}.c,Zcell{iSet}.G(:, 1:minNrOfGens)];
    end

    % obtain Zadd
    Zadd{1} = Z1.G(:, minNrOfGens+1 : end);
    for iSet = 1:length(Zcell)
        Zadd{iSet+1} = Zcell{iSet}.G(:, minNrOfGens+1 : end);
    end

    % as center is prepended
    minNrOfGens = minNrOfGens + 1;

    % compute vertex sets for each set of generators
    for iGen = 1:minNrOfGens
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
        Zmat(:,end+1 : end+n+1) = [Z_encl.c,Z_encl.G];
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

function Z = aux_unionParallelotope(Zcell)
% enclose the union of the zonotopes with a parallelotope

    % obtain matrix of points from generator matrices
    n = dim(Zcell{1});
    V = zeros(n,0);
    for i = 1:length(Zcell)
        G = generators(Zcell{i});
        V = [V, G, -G];
    end

    % compute the arithmetic mean of the vertices
    meanV = sum(V,2) / n;

    % obtain sampling matrix
    sampleMatrix = V - meanV;

    % compute the covariance matrix
    C = cov(sampleMatrix');

    % singular value decomposition
    [U,~,~] = svd(C);

    % enclose zonotopes with intervals in the transformed space
    I = interval.empty(n);
    for i = 1:length(Zcell)
        I = interval(U.' * Zcell{i}) | I;
    end

    % transform back to original space
    Z = U*zonotope(I);
end

% ------------------------------ END OF CODE ------------------------------
