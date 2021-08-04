function res = conZonotope(obj,varargin)
% conZonotope - Computes a constrained zonotope that over-approximates the 
%               conPolyZonotope object
%
% Syntax:  
%    res = conZonotope(obj)
%    res = conZonotope(obj,type)
%    res = conZonotope(obj,type,method)
%
% Inputs:
%    obj - conPolyZonotope object
%    type - algorithm used to compute conZonotope enclosure ('linearize',
%           'extend', or 'all')
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    expMat_ = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    expMat = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    cZ = conZonotope(cPZ);
%
%    figure; hold on
%    plot(cPZ,[1,2],'r','Splits',20,'Filled',true,'EdgeColor','none');
%    plot(cZ,[1,2],'b');
%    plot(zonotope(cPZ),[1,2],'g');
%
% References:
%    [1] L. Jaulin and et al. "Applied Interval Analysis", 2006 
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, interval, polyZonotope

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'all';
    method = 'linearize';
    
    if nargin > 1 && ~isempty(varargin{1})
        type = varargin{1};
    end
    if nargin > 2 && ~isempty(varargin{2})
        method = varargin{2};
    end

    % rescale the domain for the factors
    if ~strcmp(method,'none')
        obj = rescale(obj,method);
    end
    
    % check if constraints are present
    if ~isempty(obj.A)

        % enclose set with a conZonotope using the selected algorithm
        if strcmp(type,'extend')
            res = conZonotopeExtend(obj);
        elseif strcmp(type,'linearize')
            res = conZonotopeLinearize(obj);
        elseif strcmp(type,'all')
            res = conZonotopeAll(obj);
        else
           error('Wrong value for input argument "type"!'); 
        end
    
    else
        res = conZonotope(zonotope(obj));
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = conZonotopeExtend(obj)
% enclose set with a constrained zonotope based on an extension to a higher
% dimensional space

    % transform to equivalent higher-dimensional polynomial zonotope
    c = [obj.c; -obj.b];
    G = blkdiag(obj.G,obj.A);
    expMat = [obj.expMat,obj.expMat_];
    
    pZ = polyZonotope(c,G,[],expMat);
    
    % enclose higher dimensional polynomial zonotope with a zonotope
    zono = zonotope(pZ);

    % extract linear generators and constraints from zonotope enclosure
    n = dim(obj);
    G_ = generators(zono); c_ = center(zono);
    
    G = G_(1:n,:); c = c_(1:n);
    A = G_(n+1:end,:); b = -c_(n+1:end);
    
    res = conZonotope([c,G],A,b);
    
    % add independent generators
    if ~isempty(obj.Grest)
       res = res + zonotope(zeros(n,1),obj.Grest); 
    end
end

function res = conZonotopeLinearize(obj)
% enclose set with a constrained zonotope based on an enclosure of the
% nonlinear constraints with parallel hyperplanes according to Sec. 4.3.4
% in [1]

    % over-approximate the unconstrained poly. zonotope with a zonotope
    G = obj.G;
    
    temp = prod(ones(size(obj.expMat))-mod(obj.expMat,2),1);
    ind_ = find(temp == 1);
    Gquad = G(:,ind_);

    ind = find(sum(obj.expMat,1) == 1);
    G1 = zeros(size(G,1),length(obj.id));
    
    for i = 1:length(obj.id)
        index = find(obj.expMat(i,ind) == 1);
        if ~isempty(index)
            G1(:,i) = G(:,ind(index));   
        end
    end
    
    G(:,ind_) = [];
    
    c = obj.c + 0.5 * sum(Gquad,2);
    G = [G1, 0.5*Gquad G];
    
    % enclose each nonlinear constraint f(x) by two parallel hyperplanes 
    % defined as J*x + l <= f(x) <= J*x + b (see Eq. (4.93) in [1])
    jacHan = funHanJacobian(obj.A,obj.expMat_);
    funHan = @(x) funCon(x,obj.b,obj.A,obj.expMat_);
    
    temp = ones(length(obj.id),1);
    dom = interval(-temp,temp);
    
    cen = center(dom);
    J = jacHan(center(dom));
    Jint = jacHan(dom);
    
    int = funHan(cen) - J*cen + (Jint - J)*(dom - cen);   
    
    % construct the resulting constrained zonotope
    A = [J,zeros(size(J,1),size(G,2)-size(G1,2)),diag(rad(int))];
    b = -center(int);
    
    G = [G,zeros(size(G,1),size(J,1))];
    
    res = conZonotope([c,G],A,b);    
    
    % add independent generators
    if ~isempty(obj.Grest)
       res = res + zonotope(zeros(size(c)),obj.Grest); 
    end
end

function res = conZonotopeAll(obj)
% enclose set using a combination of the algorithms "extend" and
% "linearize"

    % compute enclosure with method "extend"
    res1 = conZonotopeExtend(obj);
    
    % compute enclosure with method "linearize"
    res2 = conZonotopeLinearize(obj);
   
    % intersect the results
    res = res1 & res2;
end

function f = funCon(x,b,A,expMat_)
% evaluate the constrataints of a constratined polynomial zonotope

    % initialization
    f = -b;
    
    % loop over all dependent generators
    for i = 1:size(A,2)
        f = f + A(:,i)*prod(x.^expMat_(:,i),1); 
    end
end

function J = jacobianCon(x,A,expMat_)
% compute the jacobian matrix of the polynomial constraints for a
% constrained polynomial zonotope

    % initialization
    J = 0 * repmat(x(1),[size(A{1},1),length(x)]);
    
    % loop over all variables
    for i = 1:length(expMat_)
       for j = 1:size(expMat_{i},2)
          J(:,i) = J(:,i) + A{i}(:,j)*prod(x.^expMat_{i}(:,j),1); 
       end
    end
end

function jacHan = funHanJacobian(A,expMat_)
% compute a function handle for the jacobian matrix of the constraints

    % compute exponent matrix differentiazed for each variable
    Elist = cell(size(expMat_,1),1);
    Alist = cell(size(expMat_,1),1);

    for i = 1:length(Elist)
       ind = find(expMat_(i,:) > 0);
       temp = expMat_(:,ind);
       temp(i,:) = temp(i,:) - 1;
       Elist{i} = temp;
       Alist{i} = A(:,ind) * diag(expMat_(i,ind));
    end

    % function handle for jacobian matrix
    jacHan=@(x)jacobianCon(x,Alist,Elist);
end
    
%------------- END OF CODE --------------