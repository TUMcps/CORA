function cZ = or(cZ1,varargin)
% or - Computes an over-approximation for the union of a constrained
%    zonotope and other sets
%
% Syntax:
%    cZ = or(cZ1, cZ2)
%    cZ = or(cZ1, ... , cZm)
%    cZ = or(cZ1, ... , cZm, alg)
%    cZ = or(cZ1, ... , cZm, alg, order)
%
% Inputs:
%    cZ1,...,cZm - conZonotope objects
%    alg - algorithm used to compute the union ('linProg', or 'tedrake')
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    cZ - resulting conZonotope object enclosing the union
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
%    Z = [4 2 0 0;-4 1 1 0];
%    cZ3 = conZonotope(Z,[],[]);
% 
%    % compute conZonotope that encloses the union
%    res = or(cZ1,cZ2,cZ3);
% 
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'FaceColor','r');
%    plot(cZ2,[1,2],'FaceColor','b');
%    plot(cZ3,[1,2],'FaceColor','g');
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
    alg = 'linProg';
    order = [];

    % init linprog options
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end

    % distinguish case with two and case with more sets
    if nargin == 2 || (nargin > 2 && ...
                      (isempty(varargin{2}) || ischar(varargin{2})))
                  
        S = varargin{1};

        % determine conZonotope object
        [cZ1,S] = findClassArg(cZ1,S,'conZonotope');

        % different cases depending on the class of the second set
        if isa(S,'conPolyZono')
            
            cZ = S | cZ1;
        
        elseif isa(S,'conZonotope') || isa(S,'zonotope') || ...
           isa(S,'interval') || isa(S,'zonoBundle') || ...
           isa(S,'polytope') || isnumeric(S)
            
            if ~isa(S,'conZonotope')
                S = conZonotope(S); 
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
               cZ = aux_unionLinProg({cZ1,S},order,options);
            elseif strcmp(alg,'tedrake')
               cZ = aux_unionTedrake({S,S},order,options);
            else
                throw(CORAerror('CORA:wrongValue','second',"'linProg' or 'tedrake'"));
            end

        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',cZ1,S));
        end
    
    else
        
        % parse input arguments
        Zcell = {};
        counter = [];

        for i = 2:nargin
           if isa(varargin{i-1},'contSet')
              Zcell{end+1,1} = conZonotope(varargin{i-1});
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
           cZ = aux_unionLinProg([{cZ1};Zcell],order,options);
        elseif strcmp(alg,'tedrake')
           cZ = aux_unionTedrake([{cZ1};Zcell],order,options);
        else
            throw(CORAerror('CORA:wrongValue','second',"'linProg' or 'tedrake'"));
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function Z = aux_unionTedrake(Zcell,order,options)
% compute the union by solving a linear program with respect to the
% zonotope containment constraints presented in [1]

    % construct generator matrix of the enclosing zonotope
    list = cell(length(Zcell),1);
    
    for i = 1:length(Zcell)
       temp = Zcell{i};      
       if ~isempty(temp.A)
          temp = rescale(temp);
       end
       list{i} = zonotope(temp.c,temp.G);
    end
    
    Z_ = or(list{:},'iterative',order);
    G = generators(Z_);
    
    Y = G*diag(1./sqrt(sum(G.^2,1)));
    [d,ny] = size(Y);
    Hy = [eye(ny);-eye(ny)];
    
    % construct linear constraints for each zonotope
    Aeq = [];
    beq = [];
    A = [];
    A_ = [];
    
    for i = 1:length(Zcell)
       
        % obtain generator matrix and center from the current zonotope
        [X,x,Hx,hx] = AHpolytope(Zcell{i});
        
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
    
    val = linprog(f',A,b,Aeq,beq,[],[],[],options);
    
    % construct the resulting zonotope
    ub = val(1:ny);
    lb = -val(ny+1:2*ny);
    int = interval(lb,ub);
    
    c = Y*center(int);
    G = Y * diag(rad(int));
    
    Z = conZonotope([c,G],[],[]);

end

function Z = aux_unionLinProg(Zcell,order,options)
% compute an enclosing conZonotope using linear programming. As the
% constraints for the linear program we compute the upper and lower bound
% for all zonotopes that are enclosed in the normal directions of the
% halfspace representation of the enclosing zonotope

    % construct generator matrix of the final zonotope
    list = cell(length(Zcell),1);
    
    for i = 1:length(Zcell)
       temp = Zcell{i};      
       if ~isempty(temp.A)
          temp = rescale(temp);
       end
       list{i} = zonotope(temp.c,temp.G);
    end
    
    Z_ = or(list{:},'iterative',order);
    Z_ = compact_(Z_,'zeros',eps);
    G = generators(Z_);
    
    G = G*diag(1./sqrt(sum(G.^2,1)));
    
    % compute the directions of the boundary halfspaces
    [dimG,m] = size(G);
    
    Z = zonotope([zeros(dimG,1),G]);
    Z = halfspace(Z);
    
    C = Z.halfspace.H;
    
    % compute bounds for each halfspace
    d = zeros(size(C,1),1);
    
    for i = 1:size(C,1)
       
        val = -inf;
        
        % loop over all zonotopes
        for j = 1:length(Zcell)
           
            % compute bound for the current zonotope
            valTemp = supportFunc_(Zcell{j},C(i,:)','upper');
            
            % update bound
            val = max(val,valTemp);
        end
        
        d(i) = val;
    end

    % solve linear program
    f = [zeros(dimG,1);ones(m,1)];
    
    A = [];
    b = -d;
    
    for i = 1:size(C,1)
       A = [A;-[C(i,:),abs(C(i,:)*G)]];
    end
    
    A = [A;[zeros(m,dimG),-eye(m)]];
    b = [b;zeros(m,1)];
    
    x = linprog(f',A,b,[],[],[],[],[],options);
    
    % construct final zonotope
    c = x(1:dimG);
    scal = x(dimG+1:end);
    
    Z = conZonotope([c,G*diag(scal)],[],[]);
      
end

% ------------------------------ END OF CODE ------------------------------
