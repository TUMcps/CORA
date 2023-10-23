function cZ = conZonotope(cPZ,varargin)
% conZonotope - Computes a constrained zonotope that over-approximates the 
%    constrained polynomial zonotope
%
% Syntax:
%    cZ = conZonotope(cPZ)
%    cZ = conZonotope(cPZ,type)
%    cZ = conZonotope(cPZ,type,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    type - algorithm used to compute conZonotope enclosure ('linearize',
%           'extend', or 'all')
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    EC = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    E = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    cZ = conZonotope(cPZ);
%
%    figure; hold on
%    plot(cPZ,[1,2],'FaceColor','r','Splits',10);
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

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    [type,method] = setDefaultValues({'all','linearize'},varargin); 

    % check input arguments
    inputArgsCheck({{cPZ,'att','conPolyZono'};
                    {type,'str',{'all','extend','linearize'}};
                    {method,'str',{'forwardBackward','linearize',...
                        'polynomial','interval','all','none'}}});

    % rescale the domain for the factors
    if ~strcmp(method,'none')
        cPZ = rescale(cPZ,method);
    end
    
    % check if constraints are present
    if ~isempty(cPZ.A)

        % enclose set with a conZonotope using the selected algorithm
        if strcmp(type,'extend')
            cZ = aux_conZonotopeExtend(cPZ);
        elseif strcmp(type,'linearize')
            cZ = aux_conZonotopeLinearize(cPZ);
        elseif strcmp(type,'all')
            cZ = aux_conZonotopeAll(cPZ);
        end
    
    else
        cZ = conZonotope(zonotope(cPZ));
    end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_conZonotopeExtend(cPZ)
% enclose set with a constrained zonotope based on an extension to a higher
% dimensional space

    % transform to equivalent higher-dimensional polynomial zonotope
    c = [cPZ.c; -cPZ.b];
    G = blkdiag(cPZ.G,cPZ.A);
    E = [cPZ.E,cPZ.EC];
    
    pZ = polyZonotope(c,G,[],E);
    
    % enclose higher dimensional polynomial zonotope with a zonotope
    Z = zonotope(pZ);

    % extract linear generators and constraints from zonotope enclosure
    n = dim(cPZ);
    G_ = generators(Z); c_ = center(Z);
    
    G = G_(1:n,:); c = c_(1:n);
    A = G_(n+1:end,:); b = -c_(n+1:end);
    
    res = conZonotope([c,G],A,b);
    
    % add independent generators
    if ~isempty(cPZ.GI)
       res = res + zonotope(zeros(n,1),cPZ.GI); 
    end
end

function res = aux_conZonotopeLinearize(cPZ)
% enclose set with a constrained zonotope based on an enclosure of the
% nonlinear constraints with parallel hyperplanes according to Sec. 4.3.4
% in [1]

    % over-approximate the unconstrained poly. zonotope with a zonotope
    G = cPZ.G;
    
    temp = prod(ones(size(cPZ.E))-mod(cPZ.E,2),1);
    ind_ = find(temp == 1);
    Gquad = G(:,ind_);

    ind = find(sum(cPZ.E,1) == 1);
    G1 = zeros(size(G,1),length(cPZ.id));
    
    for i = 1:length(cPZ.id)
        index = find(cPZ.E(i,ind) == 1);
        if ~isempty(index)
            G1(:,i) = G(:,ind(index));   
        end
    end
    
    G(:,ind_) = [];
    
    c = cPZ.c + 0.5 * sum(Gquad,2);
    G = [G1, 0.5*Gquad G];
    
    % enclose each nonlinear constraint f(x) by two parallel hyperplanes 
    % defined as J*x + l <= f(x) <= J*x + b (see Eq. (4.93) in [1])
    jacHan = aux_funHanJacobian(cPZ.A,cPZ.EC);
    funHan = @(x) aux_funCon(x,cPZ.b,cPZ.A,cPZ.EC);
    
    temp = ones(length(cPZ.id),1);
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
    if ~isempty(cPZ.GI)
       res = res + zonotope(zeros(size(c)),cPZ.GI); 
    end
end

function res = aux_conZonotopeAll(obj)
% enclose set using a combination of the algorithms "extend" and
% "linearize"

    % compute enclosure with method "extend"
    res1 = aux_conZonotopeExtend(obj);
    
    % compute enclosure with method "linearize"
    res2 = aux_conZonotopeLinearize(obj);
   
    % intersect the results
    res = and_(res1,res2,'exact');
end

function f = aux_funCon(x,b,A,EC)
% evaluate the constrataints of a constratined polynomial zonotope

    % initialization
    f = -b;
    
    % loop over all dependent generators
    for i = 1:size(A,2)
        f = f + A(:,i)*prod(x.^EC(:,i),1); 
    end
end

function J = aux_jacobianCon(x,A,EC)
% compute the jacobian matrix of the polynomial constraints for a
% constrained polynomial zonotope

    % initialization
    J = 0 * repmat(x(1),[size(A{1},1),length(x)]);
    
    % loop over all variables
    for i = 1:length(EC)
       for j = 1:size(EC{i},2)
          J(:,i) = J(:,i) + A{i}(:,j)*prod(x.^EC{i}(:,j),1); 
       end
    end
end

function jacHan = aux_funHanJacobian(A,EC)
% compute a function handle for the jacobian matrix of the constraints

    % compute exponent matrix differentiazed for each variable
    Elist = cell(size(EC,1),1);
    Alist = cell(size(EC,1),1);

    for i = 1:length(Elist)
       ind = find(EC(i,:) > 0);
       temp = EC(:,ind);
       temp(i,:) = temp(i,:) - 1;
       Elist{i} = temp;
       Alist{i} = A(:,ind) * diag(EC(i,ind));
    end

    % function handle for jacobian matrix
    jacHan=@(x)aux_jacobianCon(x,Alist,Elist);
end
    
% ------------------------------ END OF CODE ------------------------------
