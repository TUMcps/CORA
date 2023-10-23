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
%    f = @(x) (1 + 0.1*x(2))*(x(1) + 0.05*x(1)^3) + 0.2*x(2);
%    dom = interval([-1;-1],[1;1]);
%   
%    res = contract(f,dom,'polynomial');
%
%    figure; hold on
%    xlim([-2,2]);
%    ylim([-2,2]);
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

% Authors:       Zhuoling Li, Niklas Kochdumper
% Written:       04-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    n = length(dom);

    x = sym('x',[n,1]);

    p = f(x);

    mins = infimum(dom);
    maxs = supremum(dom);
    
    % check if function is a polynomial
    try
        for i = 1:length(p)
            sym2poly(p(i));
        end
    catch ex
        if strcmp(ex.identifier,'symbolic:sym:sym2poly:errmsg2')
            throw(CORAerror('CORA:specialError',...
                ['Contractor "polynomial" is only applicable for ', ...
                   'polynomial constraints!']));
        end
    end

    % execute contraction for each variable 
    for i = 1:n

        % loop over all constraints
        for j = 1:length(p)
        
            % test if a split is needed
            if infimum(dom(i)) < 0 && supremum(dom(i)) > 0

                 % contract first splitted domain
                 tmp = dom;
                 int1 = interval(infimum(tmp(i)),0);
                 dom(i) = int1;
                 res1_ = aux_contractBox(x(i),p(j),dom,i);

                 % contract second splitted domain
                 int2 = interval(0,supremum(tmp(i)));
                 dom(i) = int2;
                 res2_ = aux_contractBox(x(i),p(j),dom,i);
                 dom = tmp;

                 % combine the results
                 res_ = interval(min(res1_),max(res2_));

                 mins(i) = infimum(res_);
                 maxs(i) = supremum(res_);

                 if infimum(res_) < infimum(dom(i))
                     mins(i) = infimum(dom(i));
                 end

                 if infimum(res_) > supremum(dom(i))
                     maxs(i) = supremum(dom(i));
                 end
                 
                 dom = interval(mins,maxs) & dom;

                 continue
            end

            % contract overall domain
            res_ = aux_contractBox(x(i),p(j),dom,i);
            mins(i) = min(res_);
            maxs(i) = max(res_);
            dom = interval(mins,maxs) & dom;
        end
     end

    res = dom;
end


% Auxiliary functions -----------------------------------------------------

function res_ = aux_contractBox(var,poly,dom,index)
% contract domain based on extremal functions
    
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
        if (mod(d,2) == 1 && mod(k,2) == 0 && negInt == 1) || ...
           (mod(d,2) == 0 && mod(k,2) == 1 && negInt == 1)
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

    % initialize l = left bound;
    l = infimum(dom(index));

    if polyval(infs,infimum(dom(index))) <= 0 && ...
       polyval(sups,infimum(dom(index))) >= 0
        l = infimum(dom(index));
    end

    % inf(g[Y])(inf(x)) > 0
    if polyval(infs,infimum(dom(index))) > 0

        % select inf(g[Y]) for contraction
        if d < 4 
            % use explicit analytical expressions 
            l_ = roots(infs);
            if ~isempty(l_) && isreal(l_) && min(l_) > l
                l = l_;
            end
        else
            % use forward-backward contractor to find zero crossings
            int = contractForwardBackward(infHan,dom(index));

            if ~isempty(int)
                 l = infimum(int);
            end
        end
    end

    % sup(g[Y])(inf(x)) < 0
    if polyval(sups,infimum(dom(index))) < 0

        % select sup(g[Y]) for contraction
        if d < 4
            l_ = roots(sups);
            if ~isempty(l_) && isreal(l_) && min(l_) > l
                l = l_;
            end
        else
            % use forward-backward contractor to find zero crossings
            int = contractForwardBackward(supHan,dom(index));
            
            if ~isempty(int)
               l = infimum(int); 
            end
        end    
    end

    % determine the right bound of new interval

    % initialize r = right bound;
    r = supremum(dom(index));

    if polyval(infs,supremum(dom(index))) <= 0 && ...
       polyval(sups,supremum(dom(index))) >= 0
        r = supremum(dom(index));

    % inf(g[Y])(sup(x)) > 0
    elseif polyval(infs,supremum(dom(index))) > 0

        % select inf(g[Y]) for contraction
        if d < 4 
            r_ = roots(infs);
            if ~isempty(r_) && isreal(r_) && max(r_) < r
                r = max(roots(infs));
            end
        else
            % use forward-backward contractor to find zero crossings
            int = contractForwardBackward(infHan,dom(index));
            
            if ~isempty(int)
               r = supremum(int); 
            end
        end 

    % sup(g[Y])(sup(x)) < 0
    elseif polyval(sups,supremum(dom(index))) < 0

        % select sup(g[Y]) for contraction
        if d < 4
            r_ = roots(sups);
            if ~isempty(r_) && isreal(r_) && max(r_) < r
                r = max(r_);
            end
        else
            % use forward-backward contractor to find zero crossings
            int = contractForwardBackward(supHan,dom(index));
            
            if ~isempty(int)
               r = supremum(int); 
            end
        end    
    end

    res_ = [l;r];
end

% ------------------------------ END OF CODE ------------------------------
