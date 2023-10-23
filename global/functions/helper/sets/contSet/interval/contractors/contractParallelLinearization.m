function res = contractParallelLinearization(f,dom,jacHan,varargin)
% contractParallelLinearization - implementation of the parallel
%                                 linearization contractor acc. to 
%                                 Sec. 4.3.4 in [1]
%
% Syntax:
%    res = contractParallelLinearization(f,dom,jacHan)
%    res = contractParallelLinearization(f,dom,jacHan,method)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
%    jacHan - function handle for the constraint jacobian matrix
%    method - range bounding method ('interval' or 'taylm')
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) [x(1)^2 - 4*x(2); 
%              x(2)^2 - 2*x(1) + 4*x(2)];
%    dom = interval([-0.1;-0.1],[0.1;0.1]);
%   
%    res = contract(f,dom,'linearize');
%
%    figure
%    hold on
%    xlim([-0.15,0.15]);
%    ylim([-0.15,0.15]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    f = f([x1;x2]);
%    ls1 = levelSet(f(1),[x1;x2],'==');
%    ls2 = levelSet(f(2),[x1;x2],'==');
%    plot(ls1,[1,2],'b');
%    plot(ls2,[1,2],'c');
%
% References:
%    [1] L. Jaulin et al. "Applied Interval Analysis", 2006
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contract

% Authors:       Zhuoling Li, Niklas Kochdumper
% Written:       04-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    method = 'taylm';
    if nargin >= 4 && ~isempty(varargin{1})
       method = varargin{1}; 
    end

    % apply mean value theorem to enclose constraints by parallel planes
    mi = center(dom);
    A = jacHan(mi);

    if strcmp(method,'taylm')
       try
           tay = taylm(dom);
           J = interval(jacHan(tay));
       catch
           J = jacHan(dom); 
       end
    else
        J = jacHan(dom);
    end
    
    b = f(mi) - A * mi + (J - A) * (dom - mi);

    % solve linear program to compute new bounds for each variable
    A_ = [-A;A];
    b_ = [supremum(b);-infimum(b)];
    
    infi = infimum(dom);
    sup = supremum(dom);
    
    options = optimoptions('linprog','Display','off', ...
                           'ConstraintTolerance',1e-9);
    
    n = length(mi);
    f = zeros(n,1);
    
    for i = 1:n
       
        f_ = f;
        f_(i) = 1;
        
        % compute new infimum
        [~,temp] = linprog(f_,A_,b_,[],[],infi,sup,options);
        
        if ~isempty(temp)
            infi(i) = max(infi(i),temp);
        else
            poly = polytope([A_;eye(n);-eye(n)],[b_;sup;-infi]);
            
            if representsa_(poly,'emptySet',eps)
               res = [];
               return;
            end
        end
        
        % compute new supremum
        [~,temp] = linprog(-f_,A_,b_,[],[],infi,sup,options);
        
        if ~isempty(temp)
            sup(i) = min(sup(i),-temp);
        else
            poly = polytope([A_;eye(n);-eye(n)],[b_;sup;-infi]);
            
            if representsa_(poly,'emptySet',eps)
               res = [];
               return;
            end
        end
    end
    
    res = interval(infi,sup);

end
 
% ------------------------------ END OF CODE ------------------------------
