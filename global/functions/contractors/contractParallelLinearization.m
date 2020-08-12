function res = contractParallelLinearization(f,dom)
% contractParallelLinearization - implementation of the parallel
%                                 linearization contractor acc. to 
%                                 Sec. 4.3.4 in [1]
%
% Syntax:  
%    res = contractParallelLinearization(f,dom)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
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

% Author:       Zhuoling Li, Niklas Kochdumper
% Written:      04-November-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
 
    % construct symbolic variables
    n = length(dom);
    vars = sym('x',[n,1]);
    fSym = f(vars);

    % compute jacobian matrix
    jac = jacobian(fSym,vars);
    jacHan = matlabFunction(jac,'Vars',{vars});

    % apply mean value theorem
    mi = center(dom);
    A = jacHan(mi);

    try
        tay = taylm(dom);
        b = f(mi) - A * mi + (interval(jacHan(tay)) - A) * (dom - mi);
    catch
        b = f(mi) - A * mi + (jacHan(dom) - A) * (dom - mi);
    end

    % solve linear program to compute new bounds for each variable
    A_ = [-A;A];
    b_ = [supremum(b);-infimum(b)];
    
    infi = infimum(dom);
    sup = supremum(dom);
    
    options = optimoptions('linprog','Display','off');
    
    f = zeros(n,1);
    
    for i = 1:n
       
        f_ = f;
        f_(i) = 1;
        
        % compute new infimum
        [~,temp] = linprog(f_,A_,b_,[],[],infi,sup,options);
        
        if ~isempty(temp)
            infi(i) = max(infi(i),temp);
        else
            poly = mptPolytope([A_;eye(n);-eye(n)],[b_;sup;-infi]);
            
            if isempty(poly)
               res = [];
               return;
            end
        end
        
        % compute new supremum
        [~,temp] = linprog(-f_,A_,b_,[],[],infi,sup,options);
        
        if ~isempty(temp)
            sup(i) = min(sup(i),-temp);
        else
            poly = mptPolytope([A_;eye(n);-eye(n)],[b_;sup;-infi]);
            
            if isempty(poly)
               res = [];
               return;
            end
        end
    end
    
    res = interval(infi,sup);

end
 
%------------- END OF CODE --------------