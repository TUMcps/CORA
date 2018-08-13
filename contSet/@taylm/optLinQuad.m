function int = optLinQuad(obj)
% optLinQuad - bounding optimization using Linear Dominated Bounder (LDB) 
%              and Quadratic Fast Bounder (QFB)
%
% Syntax:  
%    res = optLinQuad( obj )
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
% See also: taylm, taylm/interval, optBnb, optBnbAdv
%
% References: 
%   [1] K. Makino et al. "Verified Global Optimization with Taylor Model 
%       based Range Bounders"
%   [2] K. Makino "Rigorous analysis of nonlinear motion in particle
%       accelerators"

% Author:       Niklas Kochdumper
% Written:      14-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % Implementation of the Verified Globel Optimizer concept from
    % reference paper [1] using the Linear Dominated Bounder (LDB) and the 
    % Quadratic Fast Bounder (QFB)

    % initialize all variables
    unitVec = ones(length(obj.names_of_var),1);
    dom = interval(-unitVec,unitVec);
    
    % calculate minimum
    minVal = globalMinimizer(obj,dom);
    
    % caclulate maximum
    maxVal = globalMinimizer(-obj,dom);
    
    % construct overall interval
    int = interval(minVal,-maxVal);

end


% Auxiliary Functions -----------------------------------------------------

function minVal = globalMinimizer(obj,dom)

    % initialize all variables
    domMin{1} = dom;
    tayMin{1} = obj;
    
    minVal = inf;       % global minimum   
    cutOffMin = inf;    % global minimum is guaranteed to be smaller than this value
    
    % loop over all subdomains
    while ~isempty(domMin)
        
        dom = domMin{1};
        t = tayMin{1};
            
        % split the talyor model into a part that contains the constant and
        % linear monomials, and one part the contains the high order monomials
        [objLin,objHo] = splitLinHighOrder(t);

        % compute the boundaries for both parts
        intHo = interval(objHo,'int');
        intLin = interval(objLin,'int');
        
        % update global cut-off value
        cutOffMin = min(cutOffMin,infimum(intLin)+supremum(intHo));
        
        % check for convergence and update the global minimum
        if (rad(intHo)-rad(obj.remainder)) < obj.eps/2
            minVal = min(minVal,infimum(intHo)+infimum(intLin));
        end
        
        % check if the subdomain can be eliminated
        if infimum(intLin) + infimum(intHo) > cutOffMin || (rad(intHo)-rad(obj.remainder)) < obj.eps/2
            
            domMin = domMin(2:end);
            tayMin = tayMin(2:end);            
        else

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
            if p == 0 && l > cutOffMin
                
                domMin = domMin(2:end);
                tayMin = tayMin(2:end);
            else
            
                % try to reduce the initial search domain using the LDB
                [Dmin,xMin] = calcRedDomain(objLin,intHo,'min',dom); 

                % bisect the domain along the largest dimension
                [dom1,dom2,t1,t2] = halveDomain(obj,Dmin);

                % add the subdomain that contains the minimum of the linear
                % part at position 1 in the priority queue
                if isIntersecting(dom1,xMin)
                   domMin{1} = dom1;
                   tayMin{1} = t1;
                   domMin{end+1} = dom2;
                   tayMin{end+1} = t2;
                else
                   domMin{1} = dom2;
                   tayMin{1} = t2;
                   domMin{end+1} = dom1;
                   tayMin{end+1} = t1;
                end  
            end
        end
    end
end

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

function [D,x] = calcRedDomain(objLin,intHo,type,Dorig)

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
    
    % map the reduced search domain back to the original domain
    i = mid(Dorig) + rad(Dorig) .*infimum(Dunit);
    s = mid(Dorig) + rad(Dorig) .*supremum(Dunit);
    D = interval(i,s);
    
    % map the extreme point back to the original domain
    x = mid(Dorig) + rad(Dorig) .* b1;
    
end

function [dom1,dom2,t1,t2] = halveDomain(obj,dom)

   % determine coordinate with the larges radius
   radius = rad(dom);
   [~,indCoord] = max(radius);

   % split the subdomain along the coordinate with largest radius
   infi = infimum(dom);
   sup = supremum(dom);

   temp = sup;
   temp(indCoord) = infi(indCoord) + radius(indCoord);
   dom1 = interval(infi,temp);

   temp = infi;
   temp(indCoord) = sup(indCoord) - radius(indCoord);
   dom2 = interval(temp,sup);

   % re-expand the taylor models at the center of the halved
   % domains
   t1 = reexpand(obj,dom1);
   t2 = reexpand(obj,dom2);

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

function x = gradientDescend(obj,dom)

    % determine all 2nd-order monomials
    ind = find(obj.monomials(:,1) == 2);
    mon = obj.monomials(ind,2:end);
    coeff = obj.coefficients(ind);
    
    % compute the derivative
    grad = cell(size(mon,2),1);
    
    for j = 1:size(mon,2)
        
        monTemp = mon;
        cTemp = coeff;
        
        for i = 1:size(mon,1)
            if monTemp(i,j) == 2
               cTemp(i) = 2*cTemp(i); 
            end
            monTemp(i,j) = monTemp(i,j) - 1;
        end
        
        ind = find(min(monTemp,[],2) >= 0);
        grad{j}.mon = monTemp(ind,:);
        grad{j}.coeff = cTemp(ind);
    end
    
    % constraint gradient descend
    x = mid(dom);
    dx = 0.1;
    
    while true
        
        % standard gradient descend
        G = evalGrad(grad,x);
        x_ = x - G*dx;
        
        % project the point back onto the domain
        x_ = min(max(infimum(dom),x_),x_);
        
        % check for convergence
        if norm(x-x_) < 1e-6
           x = x_;
           break; 
        end
    end
end

function G = evalGrad(grad,x)

    G = zeros(size(x));
    for i = 1:length(grad)
        G(i) = sum(prod(repmat(x,1,length(grad{i}.coeff)).^(grad{i}.mon'),1) * grad{i}.coeff);
    end
end



% Functions for debugging -------------------------------------------------

function displayRedDom2D(slopes,intHo,xMin,xMax)

    [X1,X2] = meshgrid(-1:0.1:1,-1:0.1:1);
    Z = zeros(size(X1));
    Zo = zeros(size(X1));
    Zl = zeros(size(X1));
    
    
    bound = slopes' * xMin + supremum(intHo);
    Zb = ones(size(X2)) * bound;
    
    for i = 1:size(X1,1)
        for j = 1:size(X1,2)
            x = [X1(i,j);X2(i,j)];
            Z(i,j) = slopes' * x;
            Zo(i,j) = Z(i,j) + supremum(intHo);
            Zl(i,j) = Z(i,j) + infimum(intHo);
        end
    end
    
    hold on
    s1 = surf(X1,X2,Z);
    set(s1,'EdgeColor','none');
    set(s1,'FaceColor','r');
    set(s1,'FaceAlpha',0.5);
    
    s2 = surf(X1,X2,Zo);
    set(s2,'EdgeColor','none');
    set(s2,'FaceColor','b');
    set(s2,'FaceAlpha',0.5);
    
    s3 = surf(X1,X2,Zl);
    set(s3,'EdgeColor','none');
    set(s3,'FaceColor','b');
    set(s3,'FaceAlpha',0.5);
    
    s4 = surf(X1,X2,Zb);
    set(s4,'EdgeColor','none');
    set(s4,'FaceColor','g');
    set(s4,'FaceAlpha',0.5);
    
    plot(interval(min(xMax,xMin),max(xMax,xMin)));

end

function displayGradientDescend2D(mon,coeff,xMin)

    [X1,X2] = meshgrid(-1:0.1:1,-1:0.1:1);
    Z = zeros(size(X1));
    
    for i = 1:size(X1,1)
        for j = 1:size(X1,2)
            x = [X1(i,j);X2(i,j)];
            Z(i,j) = sum(prod(repmat(x,1,length(coeff)).^(mon'),1) * coeff);
        end
    end
    
    fopt = sum(prod(repmat(xMin,1,length(coeff)).^(mon'),1) * coeff);
    
    surf(X1,X2,Z);
    hold on
    plot3(xMin(1),xMin(2),fopt,'.r','MarkerSize',20);

end

%------------ END OF CODE ------------ 