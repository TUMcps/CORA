function res = in(obj1,obj2,varargin)
% in - determines if obj2 is enclosed by the constrained zonotope obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%    res = in(obj1,obj2,type)
%    res = in(obj1,obj2,type,tol)
%
% Inputs:
%    obj1 - conZonotope object
%    obj2 - contSet object or single point
%    type - type of containment check ('exact' or 'approx')
%    tol - tolerance for the containment check
%
% Outputs:
%    res - flag specifying if set is enclosed (0 or 1)
%
% Example: 
%    % generate constrained zonotopes
%    Z = [0 2 -2 1;0 1.5 1 -1.5];
%    A = [1 1 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
% 
%    Z = [0 2 0 0;0 1 1 0];
%    A = [1 1 -1];
%    b = 0;
%    cZono2 = conZonotope(Z,A,b);
%
%    Z = [1 2 0 0;1 1 1 0];
%    A = [1 1 -1];
%    b = 0;
%    cZono3 = conZonotope(Z,A,b);
%
%    % check for containment
%    in(cZono1,cZono2)
%    in(cZono1,cZono3)
%
%    % visualization
%    figure
%    hold on
%    plot(cZono1,[1,2],'r');
%    plot(cZono2,[1,2],'b');
%
%    figure
%    hold on
%    plot(cZono1,[1,2],'r');
%    plot(cZono3,[1,2],'b');
%
% References:
%    [1] Sadraddini et. al: Linear Encodings for Polytope Containment
%        Problems, CDC 2019
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/in

% Author:       Niklas Kochdumper
% Written:      14-November-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    % parse input arguments
    type = 'exact';
    
    if nargin >= 3 && strcmp(varargin{1},'approx')
        type = 'approx';
    end
    if nargin >= 4 && ~isempty(varargin{2})
        tol = varargin{2}; n = dim(obj1);
        obj1 = obj1 + zonotope(zeros(n,1),tol*eye(n));
    end
        
    % point or point cloud in constrained zonotope containment
    if isnumeric(obj2)
        
        for i = 1:size(obj2,2)
            res = containsPoint(obj1,obj2(:,i)); 
            if res ~= 1
               return; 
            end
        end
        
    % capsule/ellipsoid in constrained zonotope containment
    elseif isa(obj2,'capsule') || isa(obj2,'ellipsoid')
        
        poly = mptPolytope(obj1);
        res = in(poly,obj2); 

    else
        
        % use the fast but over-approximative or the exact but possibly
        % slow containment check
        if strcmp(type,'exact')

            if isa(obj2,'taylm') || isa(obj2,'polyZonotope')
                error('Exact containment check not possible for this set representation!'); 
            elseif isa(obj2,'interval')
                res = in(obj1,vertices(obj2));
            else
                poly = mptPolytope(obj1);
                res = in(poly,obj2); 
            end
            
        else
            
            if isa(obj2,'taylm') || isa(obj2,'polyZonotope')
                poly = mptPolytope(obj1);
                res = in(poly,obj2); 
            else
                obj2 = conZonotope(obj2);
                res = containsSet(obj1,obj2);
            end
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = containsPoint(cZ,obj)
% use linear programming to check if a point is located inside a
% constrained zonotope

    % get object properties
    n = size(cZ.Z,1);
    p = size(cZ.Z,2)-1;
    
    c = cZ.Z(:,1);
    G = cZ.Z(:,2:end);

    % check input arguments
    if size(obj,1) ~= n || size(obj,2) ~= 1
       error('Input argument ''obj'' has the wrong format!'); 
    end

    % construct inequality constraints
    A = [eye(p);-eye(p)];
    b = ones(2*p,1);
    
    % construct equality constraints
    Aeq = [cZ.A;G];
    beq = [cZ.b;obj-c];
    
    % add slack variables
    m = size(A,1);
    meq = size(Aeq,1);
    
    A = [A,-eye(m); zeros(m,p),-eye(m)];
    b = [b;zeros(m,1)];
    Aeq = [Aeq,zeros(meq,m)];
    
    % solve linear program
    f = [zeros(p,1);ones(m,1)];
    
    options = optimoptions('linprog','display','off',...
                           'ConstraintTolerance',1e-9);
    
    [val,~,exitflag] = linprog(f',A,b,Aeq,beq,[],[],options); 
    
    % check for containment
    res = 1;
    
    if exitflag < 0 || any(val(p+1:end) > eps)
        res = 0;
    end
end

function res = containsSet(cZ1,cZ2)
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
    
    options = optimoptions('linprog','display','off');
    
    [val,~,exitflag] = linprog(f',A,b,Aeq,beq,[],[],options); 
    
    % check for containment
    res = 1;
    
    if exitflag < 0 || any(val(nx*ny+qy*qx+ny+1:end) > eps)
        res = 0;
    end
end

%------------- END OF CODE --------------