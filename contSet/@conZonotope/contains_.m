function res = contains_(cZ,S,type,tol,varargin)
% contains_ - determines if a constrained zonotope contains a set or a
%    point
%
% Syntax:
%    res = contains_(cZ,S)
%    res = contains_(cZ,S,type)
%    res = contains_(cZ,S,type,tol)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object or single point
%    type - type of containment check ('exact' or 'approx')
%    tol - tolerance for the containment check
%
% Outputs:
%    res - true/false
%
% Example: 
%    % generate constrained zonotopes
%    Z = [0 2 -2 1;0 1.5 1 -1.5];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
% 
%    Z = [0 2 0 0;0 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    Z = [1 2 0 0;1 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ3 = conZonotope(Z,A,b);
%
%    % check for containment
%    contains(cZ1,cZ2)
%    contains(cZ1,cZ3)
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ2,[1,2],'b');
%
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ3,[1,2],'b');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%    [2] JK Scott, DM Raimondo, GR Marseglia, RD Braatz: Constrained 
%        zonotopes: A new tool for set-based estimation and fault
%        detection, Automatica 69, 126-136
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper
% Written:       14-November-2019 
% Last update:   15-November-2022 (MW, return logical array for points)
%                25-November-2022 (MW, rename 'contains')
%                09-March-2023 (MA, added reference for point enclosure)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

    % enlarge constrained zonotope by tolerance
    n = dim(cZ);
    cZ = cZ + zonotope(zeros(n,1),tol*eye(n));
        
    % point or point cloud in constrained zonotope containment
    if isnumeric(S)
        
        res = false(1,size(S,2));
        for i = 1:size(S,2)
            res(i) = aux_containsPoint(cZ,S(:,i)); 
        end
        
    % capsule/ellipsoid in constrained zonotope containment
    elseif isa(S,'capsule') || isa(S,'ellipsoid')
        
        P = polytope(cZ);
        res = contains_(P,S,type,tol);

    else
        
        % use the fast but over-approximative or the exact but possibly
        % slow containment check
        if strcmp(type,'exact')

            if isa(S,'taylm') || isa(S,'polyZonotope')
                throw(CORAerror('CORA:noExactAlg',cZ,S));
            elseif isa(S,'interval')
                res = contains_(cZ,vertices(S),type,tol);
            else
                P = polytope(cZ);
                res = contains_(P,S,type,tol); 
            end
            
        else
            
            if isa(S,'taylm') || isa(S,'polyZonotope')
                P = polytope(cZ);
                res = contains_(P,S,type,tol); 
            else
                S = conZonotope(S);
                res = aux_containsSet(cZ,S);
            end
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_containsPoint(cZ,p)
% use linear programming to check if a point is located inside a
% constrained zonotope; see (20) in [2]

    % get object properties
    nrGens = size(cZ.G,2);
    
    c = cZ.c;
    G = cZ.G;

    % construct inequality constraints
    A = [eye(nrGens);-eye(nrGens)];
    b = ones(2*nrGens,1);
    
    % construct equality constraints
    Aeq = [cZ.A;G];
    beq = [cZ.b;p-c];
    
    % add slack variables
    m = size(A,1);
    meq = size(Aeq,1);
    
    A = [A,-eye(m); zeros(m,nrGens),-eye(m)];
    b = [b;zeros(m,1)];
    Aeq = [Aeq,zeros(meq,m)];
    
    % solve linear program
    f = [zeros(nrGens,1);ones(m,1)];
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off',...
                            'ConstraintTolerance',1e-9);
    end
    
    [val,~,exitflag] = linprog(f',A,b,Aeq,beq,[],[],options); 
    
    % check for containment
    res = true;
    
    if exitflag < 0 || any(val(nrGens+1:end) > eps)
        res = false;
    end
end

function res = aux_containsSet(cZ1,cZ2)
% check polytope in polytope containment according to Theorem 1 in [1]

    % convert to AH polytopes
    [Y,y,Hy,hy] = AHpolytope(cZ1);
    [X,x,Hx,hx] = AHpolytope(cZ2);
    
    nx = size(X,2);
    ny = size(Y,2);
    
    qx = length(hx);
    qy = length(hy);
    
    d = length(x);
    
    % construct constraint X = Y * T
    temp = repmat({Y},[1,nx]);
    A1 = [blkdiag(temp{:}),zeros(d*nx,qy*qx),zeros(d*nx,ny)];
    b1 = reshape(X,[d*nx,1]);

    % construct constraint y-x = Y * beta
    A2 = [zeros(d,nx*ny),zeros(d,qx*qy),Y];
    b2 = y-x;

    % construct constraint lambda * Hx = Hy * T
    Atemp = [];
    for j = 1:nx
        h = Hx(:,j);
        temp = repmat({h'},[1,qy]);
        Atemp = [Atemp;blkdiag(temp{:})];
    end

    temp = repmat({Hy},[1,nx]);
    A3 = [blkdiag(temp{:}),-Atemp,zeros(size(Atemp,1),ny)];
    b3 = zeros(size(A3,1),1);

    % construct overall equality constraints
    Aeq = [A1;A2;A3];
    beq = [b1;b2;b3];

    % construct constraint lambda * hx <= hy + Hy beta
    temp = repmat({hx'},[1,qy]);
    A1 = [zeros(qy,ny*nx),blkdiag(temp{:}),-Hy];
    b1 = hy;

    % construct constraint lambda >= 0
    A2 = [zeros(qx*qy,nx*ny),-eye(qx*qy),zeros(qx*qy,ny)];
    b2 = zeros(qx*qy,1);
    
    % constuct overall inequality constraints
    A = [A1;A2];
    b = [b1;b2];
    
    % add slack variables
    A = [A, [-eye(qy);zeros(qy*qx,qy)]; zeros(qy,size(A,2)) -eye(qy)];
    b = [b;zeros(qy,1)];
    
    Aeq = [Aeq,zeros(size(Aeq,1),qy)];
    
    % solve linear program
    f = [zeros(nx*ny+qy*qx+ny,1);ones(qy,1)];
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    
    [val,~,exitflag] = linprog(f',A,b,Aeq,beq,[],[],options); 
    
    % check for containment
    res = true;
    
    if exitflag < 0 || any(val(nx*ny+qy*qx+ny+1:end) > eps)
        res = false;
    end
end

% ------------------------------ END OF CODE ------------------------------
