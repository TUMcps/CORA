function res = contractInterval(f,dom,jacHan,varargin)
% contractInterval - contraction based on interval arithmetic
%
% Syntax:
%    res = contractInterval(f,dom,jacHan)
%    res = contractInterval(f,dom,jacHan,method)
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
%    res = contract(f,dom,'interval');
%
%    figure; hold on
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contract

% Authors:       Niklas Kochdumper
% Written:       17-December-2020 
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

    % loop over all variables
    for i = 1:length(dom)
       
        % loop over all constraints
        for j = 1:size(A,1)
            
           if abs(A(j,i)) > 1e-10
               
              % contract interval domain based on current constraints
              a = A(j,:); 
              a(i) = 0;
              temp = -(b(j) + a*dom)/A(j,i);
              dom_ = dom(i) & temp;
              
              if ~representsa_(dom_,'emptySet',eps)
                  dom(i) = dom_;
              else
                  res = []; return; 
              end
           end
        end
    end
    
    % assign ouput arguments
    res = dom;
end
 
% ------------------------------ END OF CODE ------------------------------
