function res = contractPoly(c,G,Grest,expMat,dom,varargin)
% contractPoly - contracts a interval domain to tightly enclose polynomial
%                constraints
%
% Description:
%    Contract an interval domain for polynomial constraints defined as:
%    f(x) = c + sum{i=1}prod{k=1} x(k)^E(k,i) + sum{j=1} x(j)Grest(j,:) = 0
%
% Syntax:  
%    res = contract(c,G,Grest,expMat,dom)
%    res = contract(c,G,Grest,expMat,dom,alg)
%    res = contract(c,G,Grest,expMat,dom,alg,iter)
%    res = contract(c,G,Grest,expMat,dom,alg,iter,splits)
%
% Inputs:
%    c - constant offset of the polynomial constraint
%    G - matrix of dependent generators for the polynomial constraint
%    Grest - matrix of independent generators for the polynomial constraint
%    expMat - exponent matrix for the polynomial constraint
%    dom - initial domain (class: interval)
%    alg - algorithm used for contraction. The available algorithms are 
%          'forwardBackward' (see Sec. 4.2.4 in [1]),  'linearize' 
%          (see Sec. 4.3.4 in [1]), 'polynomial' (see [2]), 'interval', and 
%          'all' (all contractors together)
%    iter - number of iteration (integer > 0 or 'fixpoint')
%    splits - number of recursive splits (integer > 0)
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) x(1)^2 + x(2)^2 - 4;
%    dom = interval([1;1],[3;3]);
%   
%    res = contract(f,dom);
%
%    figure; hold on;
%    xlim([-3,3]); ylim([-3,3]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    ls = levelSet(f([x1;x2]),[x1;x2],'==');
%    plot(ls,[1,2],'b');   
%
% References:
%    [1] L. Jaulin et al. "Applied Interval Analysis", 2006
%    [2] G. Trombettoni et al. "A Box-Consistency Contractor Based on 
%        Extremal Functions", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contractForwardBackward, contractParallelLinearization

% Author:       Niklas Kochdumper
% Written:      04-November-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    alg = 'forwardBackward';
    iter = 1;
    splits = [];
    jacHan = [];
    if nargin >= 6 && ~isempty(varargin{1})
       alg = varargin{1}; 
    end
    if nargin >= 7 && ~isempty(varargin{2})
       iter = varargin{2}; 
    end
    if nargin >= 8 && ~isempty(varargin{3})
       splits = varargin{3}; 
    end    
    if nargin >= 9 && ~isempty(varargin{4})
       jacHan = varargin{4}; 
    end
    
    % check user input
    if ischar(iter) 
        if strcmp(iter,'fixpoint')
            iter = 10000;
        else
            error('Wrong value for input argument "iter"!');
        end
    end

    if ~ismember(alg,{'forwardBackward','linearize','polynomial', ... 
                      'interval','all'})
        error('Wrong value for input argument "alg"!');           
    end
    
    % function handle for polynomial constrained function
    f = @(x) funcPoly(x,c,G,Grest,expMat);
    
    % precompute jacobian matrix
    if isempty(jacHan) && ismember(alg,{'linearize','interval','all'})
        
        % compute exponent matrix differntiazed for each variable
        Elist = cell(size(expMat,1),1);
        Glist = cell(size(expMat,1),1);
        
        for i = 1:length(Elist)
           ind = find(expMat(i,:) > 0);
           temp = expMat(:,ind);
           temp(i,:) = temp(i,:) - 1;
           Elist{i} = temp;
           Glist{i} = G(:,ind) * diag(expMat(i,ind));
        end
        
        % function handle for jacobian matrix
        jacHan = @(x) jacobianPoly(x,Glist,Grest,Elist,size(G,1), ...
                                   size(expMat,1));
    end

    % splitting of intervals considered or not
    if isempty(splits)                                  % no splitting
    
        % iteratively contract the domain
        dom_ = dom;

        for i = 1:iter

            % contract the domain using the selected algorithm
            if strcmp(alg,'forwardBackward')
                dom = contractForwardBackward(f,dom);
            elseif strcmp(alg,'linearize')
                dom = contractParallelLinearization(f,dom,jacHan,'interval');
            elseif strcmp(alg,'polynomial')
                dom = contractPolyBoxRevisePoly(c,G,Grest,expMat,dom);
            elseif strcmp(alg,'interval')
                dom = contractInterval(f,dom,jacHan,'interval');
            elseif strcmp(alg,'all')
                dom = contractForwardBackward(f,dom);
                if ~isempty(dom)
                    dom = contractInterval(f,dom,jacHan,'interval');
                    if ~isempty(dom)
                        dom = contractParallelLinearization(f,dom,jacHan,'interval');
                        if ~isempty(dom)
                           dom = contractPolyBoxRevisePoly(c,G,Grest, ...
                                                           expMat,dom); 
                        end
                    end
                end
            end
            
            % check if set is empty
            if isempty(dom)
               res = [];
               return;
            end

            % check for convergence
            if all(abs(infimum(dom)-infimum(dom_)) < eps) && ...
               all(abs(supremum(dom)-supremum(dom_)) < eps)
                break;
            else
                dom_ = dom;
            end
        end

        res = dom;
        
    else                                                % splitting
       
        % initialization
        list = {dom};
        
        % loop over the number of recursive splits
        for i = 1:splits
            
            list_ = cell(2*length(list),1);
            counter = 1;
            
            % loop over all sets in the list
            for j = 1:length(list)
                
                % determine the best dimension to split and split the 
                % domain along the determined dimension
                domSplit = bestSplit(f,list{j});
                
                % loop over all splitted domains
                for k = 1:length(domSplit)
                    
                    % check if the domain is empty
                    temp = f(domSplit{k});
                    
                    p = zeros(length(temp),1);
                    
                    if ~in(temp,p)
                       continue; 
                    end
                
                    % contract the splitted domain
                    domTemp = contractPoly(c,G,Grest,expMat, ...
                                           domSplit{k},alg,iter,[],jacHan);
                
                    % update the queue
                    if ~isempty(domTemp)
                        list_{counter} = domTemp;
                        counter = counter + 1;  
                    end
                end
            end
            
            list = list_(1:counter-1);
        end
        
        % unite all the contracted intervals
        if ~isempty(list)
            res = list{1};

            for i = 2:length(list)
               res = res | list{i}; 
            end
        else
            res = [];
        end
    end
end


% Auxiliary Functions -----------------------------------------------------

function J = jacobianPoly(x,G,Grest,expMat,n,p)
% compute the jacobian matrix for the polynomial constraint

    % initialization
    J = 0 * repmat(x(1),[n,length(expMat)]);
    x_ = x(1:p);
    
    % loop over all variables
    for i = 1:length(expMat)
       for j = 1:size(expMat{i},2)
          J(:,i) = J(:,i) + G{i}(:,j)*prod(x_.^expMat{i}(:,j),1); 
       end
    end
    
    % consider independent generators
    J = [J,Grest];
end

function f = funcPoly(x,c,G,Grest,expMat)
% evaluate the polynomial function at the point x

    % initialization
    f = c;
    n = size(expMat,1);
    x1 = x(1:n); x2 = x(n+1:end);
    
    % loop over all dependent generators
    for i = 1:size(G,2)
        f = f + G(:,i)*prod(x1.^expMat(:,i),1); 
    end

    % add independent generators
    if ~isempty(Grest)
        f = f + Grest*x2;
    end
end

function list = bestSplit(f,dom)
% determine the dimension along which it is best to split the domain

    n = length(dom);
    vol = zeros(n,1);

    % loop over all dimensions
    for i = 1:n
       
        % split the domain at the current dimension
        dom_ = split(dom,i);
        
        % evaluate constraint function for the one of the splitted sets
        val = f(dom_{1});
        
        % evaluate the quality of the split
        vol(i) = volume(val);
    end
    
    % determine best dimension to split
    [~,ind] = min(vol);
    
    % split the domain at the determined dimension
    list = split(dom,ind);
end

function res = contractPolyBoxRevisePoly(c,G,Grest,expMat,dom)
% Implementation of the contractor in [2] that is based on extremal 
% functions, where we use a specialized implementation for polynomial
% constraints to speed-up the compuations

    % initialization
    n = length(dom);

    % execute contraction for each variable 
    for i = 1:n

        % loop over all constraints
        for j = 1:length(c)
        
            % test if a split is needed
            if infimum(dom(i)) < 0 && supremum(dom(i)) > 0

                 % contract first splitted domain
                 tmp = dom;
                 int1 = interval(infimum(tmp(i)),0);
                 dom(i) = int1;
                 if ~isempty(Grest)
                     res1_ = contractBox(c(j),G(j,:),Grest(j,:),expMat,dom,i);
                 else
                     res1_ = contractBox(c(j),G(j,:),[],expMat,dom,i);
                 end

                 % contract second splitted domain
                 int2 = interval(0,supremum(tmp(i)));
                 dom(i) = int2;
                 if ~isempty(Grest)
                     res2_ = contractBox(c(j),G(j,:),Grest(j,:),expMat,dom,i);
                 else
                     res2_ = contractBox(c(j),G(j,:),[],expMat,dom,i);
                 end
                 dom = tmp;

                 % combine the results
                 res_ = interval(min(res1_),max(res2_));
                 temp = dom(i) & res_;
                 
                 if ~isempty(temp)
                    dom(i) = temp; 
                 else
                    res = []; return;
                 end

                 continue
            end

            % contract overall domain
            if ~isempty(Grest)
                res_ = contractBox(c(j),G(j,:),Grest(j,:),expMat,dom,i);
            else
                res_ = contractBox(c(j),G(j,:),[],expMat,dom,i);
            end
            
            temp = dom(i) & interval(min(res_),max(res_));
                 
             if ~isempty(temp)
                dom(i) = temp; 
             else
                res = []; return;
             end
        end
     end

    res = dom;
end

function res_ = contractBox(c,G,Grest,expMat,dom,index)
% contract domain based on extremal functions
    
    % independent vs. dependent variable
    if index <= size(expMat,1)                  % dependent variable

        % degree of the polynomial
        d = max(expMat(index,:));

        % initialization
        infs = zeros(1,d+1);
        sups = zeros(1,d+1);
        
        infs(1) = c; sups(1) = c; 

        if supremum(dom(index)) <= 0
            negInt = 1; 
        else
            negInt = 0;
        end

        % determine extrema for each coefficient function    
        for k = 0:d

            % find polynomial expression for the current coefficient
            ind = find(expMat(index,:) == k);
            
            % if the coefficient is 0, just skip the step
            if isempty(ind)
                continue
            end

            % compute interval range for the coefficient
            E = expMat(:,ind);
            E(index,:) = zeros(1,size(E,2));
            G_ = G(:,ind);
            
            g = funcPoly(dom,zeros(size(c)),G_,[],E);

            if isa(g,'double') == 1
                infs(k+1) = g;
                sups(k+1) = g;
                continue
            end 

            % if negative value is to be plugged in and the exponent is odd
            if mod(k,2) == 1 && negInt == 1
                infs(k+1) = infs(k+1) + supremum(g);
                sups(k+1) = sups(k+1) + infimum(g);
                continue
            end
            
            infs(k+1) = infs(k+1) + infimum(g);
            sups(k+1) = sups(k+1) + supremum(g);    
        end

        infs = fliplr(infs);
        sups = fliplr(sups);
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
        
    else                                    % independent variable
        
        % get coefficent in front of the current variable
        ind = index-size(expMat,1);
        coeff = Grest(:,ind);
        
        % check if domain can be updated
        if abs(coeff) > eps
        
            % set generator for current variable to 0
            Grest(:,ind) = zeros(size(Grest,1),1);

            % compute bounds for function without current variable
            g = -funcPoly(dom,c,G,Grest,expMat)/coeff;

            % extract new upper and lower bound from value
            res_ = [infimum(g);supremum(g)];
        
        else
            res_ = [infimum(dom(index)); supremum(dom(index))];
        end
    end
end

%------------- END OF CODE --------------