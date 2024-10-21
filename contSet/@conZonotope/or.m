function S_out = or(cZ,S,varargin)
% or - Computes an over-approximation for the union of a constrained
%    zonotope and other sets
%
% Syntax:
%    S_out = cZ | S
%    S_out = or(cZ,S)
%    S_out = or(cZ,S,alg,order)
%
% Inputs:
%    cZ - conZonotope objects
%    S - contSet object, numeric, cell array
%    alg - algorithm used to compute the union ('linprog' or 'tedrake')
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    S_out - convex hull enclosing the union
%
% Example: 
%    % create constrained zonotopes
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
% 
%    Z = [4 2 0 0;4 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
% 
%    % compute conZonotope that encloses the union
%    res = or(cZ1,cZ2);
% 
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'FaceColor','r');
%    plot(cZ2,[1,2],'FaceColor','b');
%    plot(res,[1,2],'k');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or

% Authors:       Niklas Kochdumper
% Written:       14-November-2019
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[alg,order] = setDefaultValues({'linprog',[]},varargin);

% check input arguments
inputArgsCheck({{cZ,'att',{'conZonotope','numeric'}},...
                {S,'att',{'contSet','numeric'}},...
                {alg,'str',{'linprog','tedrake'}},...
                {order,'att','numeric'}});

% ensure that numeric is second input argument
[cZ,S] = reorderNumeric(cZ,S);

% check dimensions
equalDimCheck(cZ,S);

% write all into a cell-array
if iscell(S)
    S_cell = [{cZ}; S];
else
    S_cell = {cZ; S};
end

% only some class allowed
if ~all(cellfun(@(S) isa(S,'conZonotope') || isa(S,'zonotope') || ...
        isa(S,'interval') || isa(S,'zonoBundle') || isnumeric(S),...
        S_cell,'UniformOutput',true))
    throw(CORAerror('CORA:noops',Z,S));
end

% convert all to constrained zonotopes
S_cell = cellfun(@(S) conZonotope(S),S_cell,'UniformOutput',false);
    
% compute over-approximation of the union with the selected algorithm
if strcmp(alg,'linprog')
    S_out = aux_unionLinprog(S_cell,order);
elseif strcmp(alg,'tedrake')
    S_out = aux_unionTedrake(S_cell,order);
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_unionTedrake(Zcell,order)
% compute the union by solving a linear program with respect to the
% zonotope containment constraints presented in [1]

    % construct generator matrix of the enclosing zonotope
    nrZ = length(Zcell);
    list = cell(nrZ,1);
    for i = 1:nrZ
        Zcell_i = Zcell{i};      
        if ~isempty(Zcell_i.A)
            Zcell_i = rescale(Zcell_i);
        end
        list{i} = zonotope(Zcell_i.c,Zcell_i.G);
    end
    Z_ = or(list{:},'iterative',order);
    G = generators(Z_);
    
    Y = G*diag(1./sqrt(sum(G.^2,1)));
    [d,ny] = size(Y);
    Hy = [eye(ny);-eye(ny)];
    
    % construct linear constraints for each zonotope
    Aeq = []; beq = [];
    A = []; A_ = [];
    
    for i = 1:nrZ
        
        % obtain generator matrix and center from the current zonotope
        [X,x,Hx,hx] = priv_AHpolytope(Zcell{i});
        
        nx = size(X,2);
        px = size(Hx,1);
        
        % construct constraint X = Y * T
        temp = repmat({Y},[1,nx]);
        A1 = [blkdiag(temp{:}),zeros(d*nx,2*px*ny),zeros(d*nx,ny)];
        b1 = reshape(X,[d*nx,1]);

        % construct constraint x = Y * beta
        A2 = [zeros(d,nx*ny),zeros(d,2*px*ny),Y];
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
        A2 = [zeros(2*px*ny,nx*ny),-eye(2*px*ny),zeros(2*px*ny,ny)];
        A_ = [A_;zeros(2*px*ny,2*ny)];
        
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
    I = interval(-val(ny+1:2*ny),val(1:ny));
    c = Y*center(I);
    G = Y * diag(rad(I));
    Z = conZonotope(c,G);
end

function Z = aux_unionLinprog(Zcell,order)
% compute an enclosing conZonotope using linear programming. As the
% constraints for the linear program we compute the upper and lower bound
% for all zonotopes that are enclosed in the normal directions of the
% halfspace representation of the enclosing zonotope

    % construct generator matrix of the final zonotope
    nrZ = length(Zcell);
    list = cell(nrZ,1);
    for i = 1:nrZ
        Zcell_i = Zcell{i};      
        if ~isempty(Zcell_i.A)
            Zcell_i = rescale(Zcell_i);
        end
        list{i} = zonotope(Zcell_i.c,Zcell_i.G);
    end
    
    Z_ = or(list{:},'iterative',order);
    Z_ = compact_(Z_,'zeros',eps);
    G = generators(Z_);
    
    G = G*diag(1./sqrt(sum(G.^2,1)));
    [n,nrGen] = size(G);
    
    % compute the directions of the boundary halfspaces
    P = polytope(Z_ - Z_.c);
    nrIneq = size(P.A,1);

    % compute bounds for each halfspace
    val = zeros(nrIneq,nrZ);
    for i = 1:nrIneq
        % loop over all zonotopes
        for j = 1:nrZ
            % compute bound for the current zonotope
            val(i,j) = supportFunc_(Zcell{j},P.A(i,:)','upper');
        end
    end
    d = max(val,[],2);

    % solve linear program
    f = [zeros(n,1);ones(nrGen,1)];

    problem.f = f';
    problem.Aineq = [-P.A, -abs(P.A*G); zeros(nrGen,n),-eye(nrGen)];
    problem.bineq = [-d;zeros(nrGen,1)];
    problem.Aeq = [];
    problem.beq = [];
    problem.lb = [];
    problem.ub = [];
    
    x = CORAlinprog(problem);
    
    % construct final zonotope
    c = x(1:n);
    scal = x(n+1:end);
    Z = conZonotope(c,G*diag(scal));
end

% ------------------------------ END OF CODE ------------------------------
