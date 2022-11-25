function p = randPoint(obj,varargin)
% randPoint - generates a random point within a constrained polynomial 
%             zonotope
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%
% Inputs:
%    obj - conPolyZono object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point inside the constrained polynomial zonotope
%
% Example: 
%    c = [0;0];
%    G = [1 0 1;0 -2 -1];
%    expMat = [1 0 3;0 1 1;0 0 0];
%    A = [1 0.5];
%    b = 0;
%    expMat_ = [1 0;1 0;0 1];
%    Grest = [0;0.1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
%
%    points = randPoint(cPZ,10,'extreme');
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',15);
%    plot(points(1,:),points(2,:),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/randPoint

% Author:       Niklas Kochdumper
% Written:      25-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    if nargin > 1 && ~isempty(varargin{1})
       N = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
       type = varargin{2}; 
    end
    
    % generate random points
    p = zeros(dim(obj),N);
    
    if strcmp(type,'standard')
        for i = 1:N
            p(:,i) = randPointStandard(obj);
        end
    elseif strcmp(type,'extreme')
        for i = 1:N
            p(:,i) = randPointExtreme(obj);
        end
    else
        [msg,id] = errWrongInput('type');
        error(id,msg);
    end
end


% Auxiliary Functions -----------------------------------------------------

function p = randPointStandard(cPZ)
% generate a random point inside a constrained polynomial zonotope

    % try with different initial points until a suitable point is found
    for i = 1:10

        % generate random factor values
        fac = length(cPZ.id);
        a = -1 + 2*rand(fac,1);

        % compute nearest point that satisfies the constraints
        objFun = @(x) sum((x-a).^2);
        
        if ~isempty(cPZ.A)
            conFun = @(x) deal([],sum(cPZ.A.*prod(x.^cPZ.expMat_,1),2)-cPZ.b);
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
        [msg,id] = errEmptySet();
        error(id,msg); 
    end
    
    p = cPZ.c + sum(cPZ.G.*prod(a_.^cPZ.expMat,1),2);
    
    % consider the indepenent generators
    if ~isempty(cPZ.Grest)
        q = size(cPZ.Grest,2);
        b = -1 + 2*rand(q,1);

        p = p + cPZ.Grest * b;
    end
end

function p = randPointExtreme(cPZ)
% generate a random point on the boundary of a constrained polynomial
% zonotope

    Grest = cPZ.Grest;
    cPZ.Grest = [];
    res = 0;
    
    % loop until a suitable point was found
    for i = 1:10

        % draw random direction vector
        n = dim(cPZ);
        d = -1 + 2*rand(n,1);

        % find minium of dependent part along this direction
        cPZ_ = d'*cPZ;

        objFun = @(x) cPZ_.c + sum(cPZ_.G.*prod(x.^cPZ_.expMat,1),2);

        if ~isempty(cPZ.A)
            conFun = @(x) deal([],-cPZ_.b + sum(cPZ_.A.* ...
                                  prod(x.^cPZ_.expMat_,1),2));
        else
            conFun = [];
        end

        p = length(cPZ.id);
        lb = -ones(p,1); ub = ones(p,1);
        
        x0 = randPoint(interval(lb,ub));

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
        [msg,id] = errEmptySet();
        error(id,msg); 
    else
        p = cPZ.c + sum(cPZ.G.*prod(x.^cPZ.expMat,1),2);
        if isempty(p)
           [msg,id] = errEmptySet();
            error(id,msg); 
        end
    end
    
    % find maximum of the independent part along this direction
    if ~isempty(Grest)
       [~,p_] = supportFunc(zonotope(zeros(n,1),Grest),d,'lower');
       p = p + p_;
    end
end

%------------- END OF CODE --------------