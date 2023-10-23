function p = randPoint_(cPZ,N,type,varargin)
% randPoint_ - generates a random point or point cloud within a constrained
%    polynomial zonotope
%
% Syntax:
%    p = randPoint_(cPZ)
%    p = randPoint_(cPZ,N)
%    p = randPoint_(cPZ,N,type)
%
% Inputs:
%    cPZ - conPolyZono object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point inside the constrained polynomial zonotope
%
% Example: 
%    c = [0;0];
%    G = [1 0 1;0 -2 -1];
%    E = [1 0 3;0 1 1;0 0 0];
%    A = [1 0.5];
%    b = 0;
%    EC = [1 0;1 0;0 1];
%    GI = [0;0.1];
%    cPZ = conPolyZono(c,G,E,A,b,EC,GI);
%
%    points = randPoint(cPZ,10,'extreme');
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(points(1,:),points(2,:),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, polyZonotope/randPoint_

% Authors:       Niklas Kochdumper
% Written:       25-January-2020
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

    % 'all' vertices not supported
    if ischar(N) && strcmp(N,'all')
        throw(CORAerror('CORA:notSupported',...
            "Number of vertices 'all' is not supported for class conPolyZono."));
    end
        
    % initialize random points
    p = zeros(dim(cPZ),N);
    
    if strcmp(type,'standard')
        i = 1;
        while i <= N
            temp = aux_randPointStandard(cPZ);
            if ~isempty(temp)
                p(:,i) = temp;
                i = i + 1;
            end
        end
    elseif strcmp(type,'extreme')
        i = 1;
        while i <= N
            temp = aux_randPointExtreme(cPZ);
            if ~isempty(temp)
                p(:,i) = temp;
                i = i + 1;
            end
        end
    else
        throw(CORAerror('CORA:noSpecificAlg',type,cPZ));
    end
    
end


% Auxiliary functions -----------------------------------------------------

function p = aux_randPointStandard(cPZ)
% generate a random point inside a constrained polynomial zonotope

    % try with different initial points until a suitable point is found
    for i = 1:10

        % generate random factor values
        fac = length(cPZ.id);
        a = -1 + 2*rand(fac,1);

        % compute nearest point that satisfies the constraints
        objFun = @(x) sum((x-a).^2);
        
        if ~isempty(cPZ.A)
            conFun = @(x) deal([],sum(cPZ.A.*prod(x.^cPZ.EC,1),2)-cPZ.b);
        else
            conFun = []; 
        end

        lb = -ones(fac,1); ub = ones(fac,1);
        options = optimoptions('fmincon','Display','off');
    
        [a_,~,res] = fmincon(objFun,a,[],[],[],[],lb,ub,conFun,options);
        
        if res >= 1
           break; 
        end
    end
    
    if res < 1
        % set is empty
        p = []; return
    end
    
    p = cPZ.c + sum(cPZ.G.*prod(a_.^cPZ.E,1),2);
    
    % consider the indepenent generators
    if ~isempty(cPZ.GI)
        q = size(cPZ.GI,2);
        b = -1 + 2*rand(q,1);

        p = p + cPZ.GI * b;
    end
end

function p = aux_randPointExtreme(cPZ)
% generate a random point on the boundary of a constrained polynomial
% zonotope

    GI = cPZ.GI;
    cPZ.GI = [];
    res = 0;
    
    % loop until a suitable point was found
    for i = 1:10

        % draw random direction vector
        n = dim(cPZ);
        d = -1 + 2*rand(n,1);

        % find minium of dependent part along this direction
        cPZ_ = d'*cPZ;

        objFun = @(x) cPZ_.c + sum(cPZ_.G.*prod(x.^cPZ_.E,1),2);

        if ~isempty(cPZ.A)
            conFun = @(x) deal([],-cPZ_.b + sum(cPZ_.A.* ...
                                  prod(x.^cPZ_.EC,1),2));
        else
            conFun = [];
        end

        p = length(cPZ.id);
        lb = -ones(p,1); ub = ones(p,1);
        
        x0 = randPoint_(interval(lb,ub),1,'standard');

        options = optimoptions('fmincon','Display','off');
        w = warning(); warning('off');
        
        [x,~,res] = fmincon(objFun,x0,[],[],[],[],lb,ub,conFun,options);
        
        if ~res
            options = optimoptions('fmincon','Display','off', ...
                                   'Algorithm','active-set');
            [x,~,res] = fmincon(objFun,x0,[],[],[],[],lb,ub,conFun,options);             
        end
        
        warning(w);
        
        % check if suitable point could be found
        if res >= 1
            break;
        end
    end
    
    % get corresponding point in state space
    if res <= 0
        p = []; return
    else
        p = cPZ.c + sum(cPZ.G.*prod(x.^cPZ.E,1),2);
        if isempty(p)
            return
        end
    end
    
    % find maximum of the independent part along this direction
    if ~isempty(GI)
       [~,p_] = supportFunc_(zonotope(zeros(n,1),GI),d,'lower');
       p = p + p_;
    end
end

% ------------------------------ END OF CODE ------------------------------
