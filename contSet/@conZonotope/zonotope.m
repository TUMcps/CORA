function Z = zonotope(cZ,varargin)
% zonotope - over-approximates a constrained zonotope with a zonotope
%
% Syntax:
%    Z = zonotope(cZ)
%    Z = zonotope(cZ,alg)
%
% Inputs:
%    cZ - conZonotope object
%    alg - algorithm used to compute enclosure ('nullSpace' or 'reduce')
%
% Outputs:
%    res - zonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
% 
%    Z1 = zonotope(cZ,'nullSpace');
%    Z2 = zonotope(cZ,'reduce');
% 
%    figure; hold on;
%    plot(cZ,[1,2],'FaceColor',[.7 .7 .7]);
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Authors:       Niklas Kochdumper
% Written:       13-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % handle trivial cases
    if representsa_(cZ,'emptySet',1e-10)
        Z = zonotope.empty(dim(cZ));
        return;
    end

    if isempty(cZ.A)
       Z = zonotope(cZ.c,cZ.G);
       return;
    end
    
    % parse input arguments
    alg = setDefaultValues({'nullSpace'},varargin);

    % check input arguments
    inputArgsCheck({{cZ,'att','conZonotope'};
                    {alg,'str',{'nullSpace','reduce'}}});
    
    % compute over-approximation using the selected algorithm
    if strcmp(alg,'nullSpace')
        Z = aux_zonotopeNullSpace(cZ);
    elseif strcmp(alg,'reduce')
        Z = aux_zonotopeReduce(cZ);
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_zonotopeNullSpace(obj)

    % compute point satisfying all constraints with pseudo inverse
    p_ = pinv(obj.A)*obj.b;
    
    % compute null-space of constraints
    T = null(obj.A);
    
    % transform boundary constraints of the factor hypercube
    m = size(obj.A,2);
    m_ = size(T,2);
    
    A = [eye(m);-eye(m)];
    b = ones(2*m,1);
    
    A_ = A*T;
    b_ = b - A*p_;
    
    % loop over all dimensions of the transformed state space
    lb = zeros(m_,1);
    ub = zeros(m_,1);
    
    persistent options
    if isempty(options)
        options = optimoptions('linprog','display','off');
    end
    
    for i = 1:m_
        
        f = zeros(m_,1);
        f(i) = 1;
        
        % compute minimum
        [~, fval] = linprog(f',A_,b_,[],[],[],[],options);
        lb(i) = fval;
        
        % compute maximum   
        [~, fval] = linprog(-f',A_,b_,[],[],[],[],options);
        ub(i) = -fval;
    end
    
    % handle case where linprog returns very close lb/ub (single point)
    dummy = norm(ub); dummy(dummy==0) = 1;
    if norm(ub-lb)/dummy < 1e-12
        meanval = 0.5*(ub+lb);
        int = interval(meanval,meanval);
    else
        int = interval(lb,ub);
    end
    
    % compute transformation matrices
    off = p_ + T*center(int);
    S = T*diag(rad(int));
    
    % construct final zonotope
    c = obj.c + obj.G*off;
    G = obj.G*S;
    
    res = zonotope(c,G);

end

function res = aux_zonotopeReduce(obj)

    % remove all constraints of the constrained zonotope
    ng = floor((size(obj.G,2))/length(obj.c)) + 1;

    obj = reduce(obj,'girard',ng,0);

    % construct the resulting zonotope object
    res = zonotope(obj.c,obj.G);
end

% ------------------------------ END OF CODE ------------------------------
