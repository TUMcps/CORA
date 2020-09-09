function res = zonotope(obj,varargin)
% zonotope - over-approximates a constrained zonotope with a zonotope
%
% Syntax:  
%    res = zonotope(obj)
%    res = zonotope(obj,alg)
%
% Inputs:
%    obj - conZonotope object
%    alg - algorithm used to compute enclosure ('nullSpace' or 'reduce')
%
% Outputs:
%    res - zonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    zono1 = zonotope(cZono,'nullSpace');
%    zono2 = zonotope(cZono,'reduce')
%
%    hold on
%    plot(cZono,[1,2],'FaceColor',[.7 .7 .7],'Filled',true,'EdgeColor','none')
%    plot(zono1,[1,2],'b');
%    plot(zono2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  09-September-2020 (VG, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(obj)
        res = zonotope();
        return;
    end

    % parse input arguments
    alg = 'nullSpace';
    if nargin >= 2
        alg = varargin{1};
    end
    
    % compute over-approximation using the selected algorithm
    if strcmp(alg,'nullSpace')
        res = zonotopeNullSpace(obj);
    elseif strcmp(alg,'reduce')
        res = zonotopeReduce(obj);
    else
        error('Wrong value for input argument ''alg''!');
    end
end


% Auxiliary Functions -----------------------------------------------------

function res = zonotopeNullSpace(obj)

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
    
    options = optimoptions('linprog','display','off');
    
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
    
    int = interval(lb,ub);
    
    % compute transformation matrices
    off = p_ + T*center(int);
    S = T*diag(rad(int));
    
    % construct final zonotope
    c = obj.Z(:,1) + obj.Z(:,2:end)*off;
    G = obj.Z(:,2:end)*S;
    
    res = zonotope([c,G]);

end

function res = zonotopeReduce(obj)

    % remove all constraints of the constrained zonotope
    ng = floor((size(obj.Z,2)-1)/size(obj.Z,1)) + 1;

    obj = reduce(obj,'girard',ng,0);

    % construct the resulting zonotope object
    res = zonotope(obj.Z);
end

%------------- END OF CODE --------------