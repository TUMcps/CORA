function res = inverse(obj)
% inverse - Compute formula '1/obj' for Taylor models
%
% Syntax:  
%    res = inverse(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"
%   [2] M. Althoff et al. "Implementation of Taylor models in CORA 2018
%       (Tool Presentation)"
%   [3] K. Makino et al. "Verified Global Optimization with Taylor Model 
%       based Range Bounders"

% Author:       Dmitry Grebenyuk
% Written:      10-August-2017
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------

        
    if isempty(obj.monomials)    % Taylor model without polynomial part
        
        res = obj;
        res.remainder = 1/(obj.remainder);
        
    else                        % Taylor model with polynomial part
        
        if obj.monomials(1,1) == 0                 
            c_f = obj.coefficients(1);
        else
            c_f = obj.coefficients( obj.monomials(:,1) == 0 );
        end
        
        %producing T_tilda from the manual
        if isempty(c_f) || c_f == 0
            error('Function ''taylm/inverse'' is not defined for value 0!'); 
        else    
            T = obj - c_f;
            c_f = 1./c_f;
        end        
    
        % check if the input taylor model is admissible (> 0)
        rem = interval(T);
        temp = rem + 1./c_f;
        if infimum(temp) <= 0 && supremum(temp) >= 0
            error('Function ''taylm/inverse'' is not defined for value 0!'); 
        end
        
        % sum initialisation
        T1 = c_f;
        T_factor = c_f;
        
        for i = 1:obj.max_order
            T_factor = T_factor .* T .* c_f;
            T1 = T1 + (-1).^i .* T_factor;
        end

        rem = interval(T);
        
        % temporaly increase max-order to get tight bounds for interval "remPow"
        T_factor.max_order = T_factor.max_order + T.max_order;
        remPow = interval(T_factor.*T.*c_f);     % remPow = B( (T/c_f)^(k+1) )

        LagRem =  (-1).^(obj.max_order + 1) .* remPow .* ...
                   1 ./ (1 + interval(0,1) .* rem .* c_f).^(obj.max_order + 2);

        % Heuristic: if the remainder is larger than the evaluation of the
        % inverse with interval arithmetic, we assume that the remainder
        % contains a large over-approximation. In this case, we re-calculate
        % the real value of the Lagrange remainder up to a certain
        % precision with an algorithm based on splits into subdomains
        remReal = rem + 1/c_f;
        temp = 1/remReal;
        
        if infimum(LagRem) < infimum(temp) && supremum(LagRem) > supremum(temp)
            LagRem = globalBounds(1/c_f,obj.max_order,remReal,obj.eps);
        end
        
        % assemble the resulting taylor model
        res = T1;
        res.remainder = res.remainder + LagRem;
    end
end

function rem = calculateRemainder(x0,n,x)
% Equation for the exact Lagrange remainder (obtained by analytical
% integration of the Lagrange remainder term, no over-approximation
% included)

    % Implementation of equation (38) in reference paper [2]
    rem = (1./x0 - 1./ x);

    % Implementation of equation (37) in reference paper [2]
    for i = 1:n
       rem = 1/(1+i) * (x-x0).^i./x0.^(i+1) - i/(i+1) * rem; 
    end

    rem = rem * (-1)^(n+1) * (n+1);
end

function int = globalBounds(x0,n,dom,tol)
% determine the global bounds of the Lagrange remainder up to the precision 
% "tol" with interval arithmetics and splits of the original domain into
% subdomains

    % calculate minimum
    func = @(x) calculateRemainder(x0,n,x);
    minVal = globalMinimizer(func,dom,tol);
    
    % caclulate maximum
    func = @(x) -calculateRemainder(x0,n,x);
    maxVal = globalMinimizer(func,dom,tol);
    
    % construct overall interval
    int = interval(minVal,-maxVal);
end

function minVal = globalMinimizer(func,dom,tol)
% Implementation of the Verified Global Optimizer Concept from reference
% paper [3] using interval arithmetic only

    % initialize all variables
    domMin{1} = dom;
    intMin{1} = func(dom);
    
    minVal = inf;       % global minimum   
    cutOffMin = inf;    % global minimum is guaranteed to be smaller than this value
    
    % loop over all subdomains
    while ~isempty(domMin)
        
        dom = domMin{1};
        int =  intMin{1};        
        
        % update global cut-off value
        cutOffMin = min(cutOffMin,supremum(int));
        
        % check for convergence and update the global minimum
        if rad(int) < tol
            minVal = min(minVal,infimum(int));
        end
        
        % check if the subdomain can be eliminated
        if infimum(int) > cutOffMin || rad(int) < tol
            
            domMin = domMin(2:end);        
        else

            % split the interval
            dom1 = interval(infimum(dom),mid(dom));
            dom2 = interval(mid(dom),supremum(dom));
            
            % calculate new bounds
            int1 = func(dom1);
            int2 = func(dom2);

            % add the subdomain that contains the minimum of the linear
            % part at position 1 in the priority queue
            if infimum(int1) < infimum(int2)
               domMin{1} = dom1;
               intMin{1} = int1;
               domMin{end+1} = dom2;
               intMin{end+1} = int2;
            else
               domMin{1} = dom2;
               intMin{1} = int2;
               domMin{end+1} = dom1;
               intMin{end+1} = int1;
            end  
        end
   end
end


%------------ END OF CODE ------------