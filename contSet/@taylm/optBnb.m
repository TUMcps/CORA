function res = optBnb( obj )
% optBnb - a branch and bound optimization
%
% Syntax:  
%    res = optBnb( obj )
%
% Inputs:
%    obj - taylm
%
% Outputs:
%    res - interval
%
% Example:
%   f= @(x) 1 + x.^5 - x.^4;
%   x = interval(0,1);
%   t = taylm(x,10,'x');
%   T = f(t);
%
%   intReal = interval(0.91808,1)
%
%   int = interval(T)
%   intBnb = interval(T,'bnb')
%   intBnbAdv = interval(T,'bnbAdv')
%   intLinQuad = interval(T,'linQuad')
%
%   X = 0:0.01:1;
%   Y = f(X);
%   plot(X,Y);
%   xlim([0,1]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, taylm/interval, optBnbAdv, optLinQuad
%
% References: 
%   [1] M. Althoff et al. "Implementation of Taylor models in CORA 2018
%       (Tool Presentation)"

% Author:       Dmitry Grebenyuk
% Written:      23-October-2017
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % Implementation of Alogrithm 1 from reference paper [1]
    % (standard version, without reexpansion)

    tmp = length(obj.names_of_var);
    dom = interval(-ones(1, tmp), ones(1, tmp));
    
    int = 0;
    % insert the initial interval into the polynomial 
    for i = 1:length(obj.coefficients)
        exp = repmat(obj.monomials(i,2:end), size(dom, 1) ,1);
        int = int + prod(dom.^exp, 2) * obj.coefficients(i);
    end

    int_u = int;
    
    while true
               
        % find the index of the interval containing the upper bound
        [~, ind_mx] = max(supremum(int));
         % find the index of the interval containing the lower bound
        [~, ind_mn] = min(infimum(int));
        
        % halve the interval
        if ind_mx == ind_mn
            [~, ind_var] = max(rad(dom(ind_mx, :)));
            dom = halve(dom, ind_mx, ind_var);
        else
            [~, ind_var] = max(rad(dom(ind_mx, :)));
            dom = halve(dom, ind_mx, ind_var);
            
            [~, ind_var] = max(rad(dom(ind_mn, :)));
            dom = halve(dom, ind_mn, ind_var);
        end
        
        prev_int_u = int_u;
        int = 0;

        % insert the new interval into the polynomial 
        for i = 1:length(obj.coefficients)
            exp = repmat(obj.monomials(i,2:end), size(dom, 1) ,1);
            int = int + prod(dom.^exp, 2) * obj.coefficients(i);
        end
        
        int_u = unite_ints(int);
        if intcmp(prev_int_u, int_u,obj.eps)
            break;
        end
    end 
    res = int_u + obj.remainder;
    
end

function res = halve(obj, ind_bit, ind_var) % ind_bit - index of a bit ( a part of an integral)
                                            % ind_val - index of a varible
    % find a middle point
    m = mid(obj(ind_bit, ind_var));
    
    % break the chosen interval by half
    inf = infimum(obj);
    new_line_inf = inf(ind_bit, :);
    new_line_inf(:, ind_var) = m;
    
    sup = supremum(obj);
    new_line_sup = sup(ind_bit, :);
    new_line_sup(:, ind_var) = m;
    
    inf = [inf(1:ind_bit, :); new_line_inf; inf(ind_bit+1:end, :)];
    sup = [sup(1:ind_bit-1, :); new_line_sup; sup(ind_bit:end, :)];
    
    res = interval(inf, sup);
end

function res = intcmp(prev_int, int, eps)
    
    % cut off if the refinement is less than EPS
    prev_int_s = supremum(prev_int); 
    prev_int_i = infimum(prev_int);
    int_s = supremum(int);
    int_i = infimum(int);
    
    if (abs(prev_int_s - int_s) <= eps * abs(prev_int_s - prev_int_i)) && ...
          (abs(int_i - prev_int_i) <= eps * abs(prev_int_s - prev_int_i))
      res = 1;
    else
        res = 0;
    end
end

%------------- END OF CODE --------------