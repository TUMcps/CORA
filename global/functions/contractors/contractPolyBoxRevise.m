function res = contractPolyBoxRevise(f,dom)
% contractPolyBoxRevise - implementation of the contractor based on 
%                         extremal functions acc. to [1]
%
% Syntax:  
%    res = contractPolyBoxRevise(f,dom)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) x(1)^2 + x(2)^2 - 4;
%    dom = interval([1;1],[3;3]);
%   
%    res = contract(f,dom,'polyBox');
%
%    figure
%    hold on
%    xlim([-3,3]);
%    ylim([-3,3]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    ls = levelSet(f([x1;x2]),[x1;x2],'==');
%    plot(ls,[1,2],'b'); 
%
% References:
%    [1] G. Trombettoni et al. "A Box-Consistency Contractor Based on 
%        Extremal Functions", 2010
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

% number of the variables
n = length(dom);

x = sym('x',[n,1]);

p = f(x);

mins = zeros(n,1);
maxs = zeros(n,1);

% execute contraction for each variable 
for i = 1:n
        
    % test if a split is needed
     if infimum(dom(i)) < 0 && supremum(dom(i)) > 0
         tmp = dom;
         int1 = interval(infimum(tmp(i)),0);
         dom(i) = int1;
         res1_ = contractBox(x(i),p,dom,i);
         int2 = interval(0,supremum(tmp(i)));
         dom(i) = int2;
         res2_ = contractBox(x(i),p,dom,i);
         dom = tmp;
         
         res_ = interval(min(res1_),max(res2_));
                  
         mins(i) = infimum(res_);
         maxs(i) = supremum(res_);
         
         if infimum(res_) < infimum(dom(i))
             mins(i) = infimum(dom(i));
         end
         
         if infimum(res_) > supremum(dom(i))
             maxs(i) = supremum(dom(i));
         end
                     
         continue
     end
  
    res_ = contractBox(x(i),p,dom,i);
    mins(i) = min(res_);
    maxs(i) = max(res_);
 end

res = interval(mins,maxs);

end

function res_ = contractBox(var,poly,dom,index)
    
c = coeffs(poly,var,'All');

% degree of the polynomial
d = length(c);

% determine the extremal functions
infs = zeros(1,d);
sups = zeros(1,d);

n = length(dom);
x = sym('x',[n,1]);

if supremum(dom(index)) <= 0
    negInt = 1; 
else
    negInt = 0;
end
    
% determine extrema for each coefficient function    
    for k = 1:d
    
    % if the coefficient is 0, just skip the step
    if c(k) == 0
        continue
    end
        
    g = matlabFunction(c(k),'Vars',{x});
    
    if isa(g(dom),'double') == 1
        infs(k) = g(dom);
        sups(k) = g(dom);
        continue
    end 
   
    % if negative value is to be plugged in and the exponent is odd
    if (mod(d,2) == 1 && mod(k,2) == 0 && negInt == 1) || (mod(d,2) == 0 && mod(k,2) == 1 && negInt == 1)
        infs(k) = supremum(g(dom));
        sups(k) = infimum(g(dom));
        continue
    end
    
    infs(k) = infimum(g(dom));
    sups(k) = supremum(g(dom));     
	
    end

inffuc = poly2sym(infs);
infHan = matlabFunction(inffuc);
supfuc = poly2sym(sups);
supHan = matlabFunction(supfuc);

% determine the left bound of new interval

% initialize r = right bound;
l = infimum(dom(index));

if polyval(infs,infimum(dom(index))) <= 0 && polyval(sups,infimum(dom(index))) >= 0
    l = infimum(dom(index));
end

% inf(g[Y])(inf(x)) > 0
if polyval(infs,infimum(dom(index))) > 0
    
    % select inf(g[Y]) for contraction
    if d < 4 
        % use explicit analytical expressions 
        l = roots(infs);
    else
        % use HC4-revise algorithm      
        l = infimum(contractForwardBackward(infHan,dom(index)));
    end
    
    if min(l) > infimum(dom(index))
        l = min(roots(infs));
    end
end

% sup(g[Y])(inf(x)) < 0

if polyval(sups,infimum(dom(index))) < 0
    
    % select sup(g[Y]) for contraction
    if d < 4
        l_ = roots(sups);
    else
        l_ = infimum(contractForwardBackward(supHan,dom(1)));
    end    
    
    if min(l_) > l
        l = l_;
    end
end

% determine the right bound of new interval

% initialize r = right bound;
r = supremum(dom(index));

if polyval(infs,supremum(dom(index))) <= 0 && polyval(sups,supremum(dom(index))) >= 0
    r = supremum(dom(index));
    
% inf(g[Y])(sup(x)) > 0

elseif polyval(infs,supremum(dom(index))) > 0
    
    % select inf(g[Y]) for contraction
    if d < 4 
        r = roots(infs);
    else
        r = supremum(contractForwardBackward(infHan,dom(index)));
    end
    
    % if the root is greater than right bound, then no contraction needed
    if max(r) < supremum(dom(index))
        r = max(roots(infs));
    end
    
% sup(g[Y])(sup(x)) < 0
elseif polyval(sups,supremum(dom(index))) < 0
    
    % select sup(g[Y]) for contraction
    if d < 4
        r_ = roots(sups);
    else
        r_ = supremum(contractForwardBackward(supHan,dom(index)));
    end    
    
    if max(r_) < r
        r = r_;
    end
end

res_ = [l;r];

end

%------------- END OF CODE --------------