function [val,xOpt,domOpt] = globVerMinimization(func,dom,tol,varargin)
% globVerMinimization - determine the gloabl minimum of a function on a 
%                       search domain with a certain precision
%
% Syntax:  
%    [val,xOpt,domOpt] = globVerMinimization(func,dom,tol)
%    [val,xOpt,domOpt] = globVerMinimization(func,dom,tol,max_order,opt_method,eps,tolerance,objTay,objAff)
%
% Inputs:
%    func - function that is minimizes (provided as a function handle)
%    dom - multi-dimensional search domain (class: interval)
%    tol - tolerance for the determined value of the minimum. The real
%          minimum "min_real" is located inside the following interval:
%          min_real \in [val, val + tol]
%    max_order - the maximal polynomial order of a monomial stored in a 
%                polynomial part of a taylor model
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models during the calculation of the initial
%                 taylor model object
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of taylor models
%                  'linQuad': optimization with Linear Dominated Bounder (LDB)
%                             and Quadratic Fast Bounder (QFB)
%   eps - precision for the selected optimization method for the talyor 
%         over-approximation with intervals (opt_method = 'bnb', 
%         opt_method = 'bnbAdv' and opt_method = 'linQuad')
%   tolerance - taylor model monomials with coefficients smaller than this
%               value are moved to the remainder
%   objTay - taylor model object resulting from the evaluation of the
%            function "func" with taylor models (class: taylm) 
%   objAff - affine object resulting from the evaluation of the
%            function "func" with affine arithmetic (class: affine) 
%
% Outputs:
%    val - determined lower bound for the minimum of the function
%    xOpt - best guess for the point at which the function reaches it's
%           minimum
%    domOpt - domain in which the determined lower bound of the function 
%             minimum is located  
%
% Example:
%   % 2D Beale function with minimum 0 at [3,0.5]
%   f=@(x) (1.5 - x(1)*(1-x(2))).^2 + (2.25 - x(1)*(1-x(2)^2))^2 + (2.625 - x(1)*(1-x(2)^3))^2;
%   x = interval([-4.5;-4.5],[4.5;4.5]);
%   [val,xOpt,domOpt] = globVerMinimization(f,x,1e-5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: globVerBounds, optLinQuad
%
% References: 
%   [1] K. Makino et al. "Verified Global Optimization with Taylor Model 
%       based Range Bounders"
%   [2] K. Makino "Rigorous analysis of nonlinear motion in particle
%       accelerators"

% Author:       Niklas Kochdumper
% Written:      16-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % Implementation of the Verified Global Optimizer concept from
    % reference paper [1]
    

    % Initialization ------------------------------------------------------

    % default values
    max_order = 10;
    opt_method = 'int';
    eps = 0.001;
    tolerance = 1e-8;
    
    % parse input arguments
    if nargin < 3
       error('Wrong syntax!. Type "help globVerMinimization" for help.'); 
    end
    if nargin >= 4 && ~isempty(varargin{1})
       max_order = varargin{1}; 
    end
    if nargin >= 5 && ~isempty(varargin{2})
       opt_method = varargin{2}; 
    end
    if nargin >= 6 && ~isempty(varargin{3})
       eps = varargin{3}; 
    end
    if nargin >= 7 && ~isempty(varargin{4})
       tolerance = varargin{4}; 
    end
    
    % create taylor models
    if nargin >= 8 && ~isempty(varargin{5})
        T = varargin{5};
    else
        t = taylm(dom,max_order,'x',opt_method,eps,tolerance);
        T = func(t);
    end
    
    % create affine arithmetic objects
    if nargin == 9 && ~isempty(varargin{6})
        A = varargin{6};
    else
        a = affine(dom,'x',opt_method,eps,tolerance);
        A = func(a);
    end
    
    % initialize all variables
    unitVec = ones(length(dom),1);
    domTay{1} = interval(-unitVec,unitVec);
    domReal{1} = dom;
    aff{1} = A;
    tay{1} = T;
    
    minVal = inf;    % global minimum   
    cutOff = inf;    % global minimum is guaranteed to be smaller than this value
    
    
    
    
    % Main Algorithm ------------------------------------------------------
    
    % loop over all subdomains
    while ~isempty(domTay)
        
        domR = domReal{1};
        domT = domTay{1};
        a = aff{1};
        t = tay{1};
        
        elim = 0;
        
        % loop over all methods that are tested
        for i = 1:3
            
           switch i
               
               % Interval arithmetic
               case 1
               
                    temp = func(domR);
                    maxTemp = supremum(temp);
                    minTemp = infimum(temp);
                   
               % Affine arithmetic    
               case 2
                   
                    temp = interval(a,'int');
                    maxTemp = supremum(temp);
                    minTemp = infimum(temp);
                   
               % Taylor model (Linear Dominated Bounder)
               case 3
            
                    % split the talyor model into a part that contains the constant and
                    % linear monomials, and one part the contains the high order monomials
                    [objLin,objHo] = splitLinHighOrder(t);

                    % compute the boundaries for both parts
                    intHo = interval(objHo,'int');
                    intLin = interval(objLin,'int');
                    
                    % compute the lower and upper bound for the minimum
                    maxTemp = infimum(intLin) + supremum(intHo);
                    minTemp = infimum(intLin) + infimum(intHo);    
           end
           
           % update global cut-off value
           cutOff = min(cutOff,maxTemp);

           % check for convergence and update the global minimum
           if (maxTemp-minTemp) < tol/2
               if minTemp < minVal
                   minVal = minTemp;
                   domOpt = domR;
               end
           end

           % check if the subdomain can be eliminated
           if minTemp > cutOff || (maxTemp - minTemp) < tol/2

                domTay = domTay(2:end);
                domReal = domReal(2:end);
                aff = aff(2:end);
                tay = tay(2:end); 
                
                elim = 1;
                break;
           end
        end
        
                    
        if elim == 0

            % Compute a bound using the Quadratic Fast Bounder (QFB)
            H = hessian(t);
            [~,p] = chol(H);
            
            if p == 0       % Hessian matrix is positive definite
                
                % create domain of the taylor model
                unitVec = ones(length(t.names_of_var),1);
                int = interval(-unitVec,unitVec);
                
                % determine point the minimizes the quadratic part (since
                % the domain of the taylor models is always scaled to
                % [-1,1]^n, the minimum of the quadratic part is always
                % located at 0)
                x0 = zeros(length(t.names_of_var),1);
                
                % construct taylor model for term Q
                x = taylm(int,t.max_order,transpose(t.names_of_var),t.opt_method,t.eps,t.tolerance);
                Q = 0.5* transpose(x-x0) * H * (x-x0);
                
                % caclulate the lower bound
                l = infimum(interval(t-Q,'int'));           
            end
            
            % check if the subdomain can be eliminated
            if p == 0 && l > cutOff
                
                domTay = domTay(2:end);
                domReal = domReal(2:end);
                aff = aff(2:end);
                tay = tay(2:end);
                
            else
            
                % try to reduce the initial search domain using the LDB
                [Dmin,DminR,xMin] = calcRedDomain(objLin,intHo,'min',domT,domR);

                % bisect the domain along the largest dimension
                [dom1,dom2,domR1,domR2,t1,t2,a1,a2] = halveDomain(T,A,Dmin,DminR);

                % add the subdomain that contains the minimum of the linear
                % part at position 1 in the priority queue
                if isIntersecting(dom1,xMin)
                   domTay{1} = dom1;
                   domReal{1} = domR1;
                   tay{1} = t1;
                   aff{1} = a1;
                   
                   domTay{end+1} = dom2;
                   domReal{end+1} = domR2;
                   tay{end+1} = t2;
                   aff{end+1} = a2;
                else
                   domTay{1} = dom2;
                   domReal{1} = domR2;
                   tay{1} = t2;
                   aff{1} = a2;
                   
                   domTay{end+1} = dom1;
                   domReal{end+1} = domR1;
                   tay{end+1} = t1;
                   aff{end+1} = a1;
                end  
            end
        end
    end
    
    % assign the output arguments
    val = minVal;
    xOpt = mid(domOpt);
end




% Auxiliary functions -----------------------------------------------------

function [objLin,objHo] = splitLinHighOrder(obj)

    % find all monomials that are constant or linear
    ind = find(obj.monomials(:,1) <= 1);
    
    % construct the linear taylor model
    objLin = obj;
    objLin.monomials = obj.monomials(ind,:);
    objLin.coefficients = obj.coefficients(ind);
    objLin.remainder = interval(0,0);
    
    % construct the taylor model containing the higher order terms
    objHo = obj;
    objHo.monomials(ind,:) = [];
    objHo.coefficients(ind) = [];
    
end

function [D,Dr,x] = calcRedDomain(objLin,intHo,type,Dorig,DorigR)

    % extract the slopes of the linear part
    slopes = extractSlopes(objLin);

    % find the vertex of the domain that contains the minimum of the 
    % linear function
    if strcmp(type,'min')
        b1 = sign(-slopes);
    else
        b1 = sign(slopes);
    end
    ind = find(b1 == 0);
    b1(ind) = ones(length(ind),1);

    % calculate the reduced search domain for all dimensions
    b2 = -b1;

    for i = 1:length(b1)
        if slopes(i) ~= 0
           temp = (2*rad(intHo))/slopes(i) + b1(i);
           if temp > -1 && temp < 1
              b2(i) = temp; 
           end
        end 
    end

    % generate the reduced search domain (subset of unit interval)
    Dunit = interval(min(b1,b2),max(b1,b2));
    
    % map the reduced search domain back to the original domain (taylor)
    try
    i = mid(Dorig) + rad(Dorig) .*infimum(Dunit);
    catch
        test = 1;
    end
    s = mid(Dorig) + rad(Dorig) .*supremum(Dunit);
    D = interval(i,s);
    
    % map the reduced search domain back to the original domain
    i = mid(DorigR) + rad(DorigR) .*infimum(Dunit);
    s = mid(DorigR) + rad(DorigR) .*supremum(Dunit);
    Dr = interval(i,s);
    
    % map the extreme point back to the original domain
    x = mid(Dorig) + rad(Dorig) .* b1;
    
end

function [dom1,dom2,domR1,domR2,t1,t2,a1,a2] = halveDomain(T,A,dom,domR)

   % determine coordinate with the larges radius
   radius = rad(domR);
   [~,indCoord] = max(radius);

   % split the real subdomain along the coordinate with the largest radius
   infi = infimum(domR);
   sup = supremum(domR);

   temp = sup;
   temp(indCoord) = infi(indCoord) + radius(indCoord);
   domR1 = interval(infi,temp);

   temp = infi;
   temp(indCoord) = sup(indCoord) - radius(indCoord);
   domR2 = interval(temp,sup);
   
   % split the subdomain for the taylor models along the coordinate with 
   % largest radius
   infi = infimum(dom);
   sup = supremum(dom);
   radius = rad(dom);

   temp = sup;
   temp(indCoord) = infi(indCoord) + radius(indCoord);
   dom1 = interval(infi,temp);

   temp = infi;
   temp(indCoord) = sup(indCoord) - radius(indCoord);
   dom2 = interval(temp,sup);

   % adapt the taylor model domains to the possibly alternated order of 
   % variables
   [~,ind] = sort(T.names_of_var);
   dom1_ = dom1;
   dom1_(ind) = dom1;
   dom2_ = dom2;
   dom2_(ind) = dom2;
   
   % re-expand the taylor models at the center of the halved
   % domains
   t1 = reexpand(T,dom1_);
   t2 = reexpand(T,dom2_);
   
   % re-expand the affine models at the center of the halved
   % domains
   a1 = reexpand(A,dom1_);
   a2 = reexpand(A,dom2_);

end

function slopes = extractSlopes(obj)

    % extract the slopes of the linear part of the taylor model
    slopes = zeros(length(obj.names_of_var),1);
    
    for i = 1:length(slopes)
        
       ind = find(obj.monomials(:,1+i) == 1);
       
       if ~isempty(ind)
          slopes(i) = obj.coefficients(ind); 
       end   
    end
    
    % sort the slopes according to the possibly differing order of
    % variables in the taylor model
    [~,ind] = sort(obj.names_of_var);
    slopes = slopes(ind);
end

function H = hessian(obj)
    
    % determine all 2nd-order monomials
    ind = find(obj.monomials(:,1) == 2);
    mon = obj.monomials(ind,2:end);
    coeff = obj.coefficients(ind);
    
    % construct hessian matrix
    H = zeros(size(mon,2));
    for i = 1:length(coeff)
        ind = find(mon(i,:) > 0);
        if length(ind) == 1
           H(ind,ind) = 2 * coeff(i); 
        else
           H(ind(1),ind(2)) = coeff(i);
           H(ind(2),ind(1)) = coeff(i);
        end   
    end
end
    
%------------ END OF CODE ------------