function val = supportFunc(pZ,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a polynomial zonotope
%               object along a certain direction
%
% Syntax:  
%    val = supportFunc(pZ,dir)
%    val = supportFunc(pZ,dir,type)
%    val = supportFunc(pZ,dir,type,method)
%    val = supportFunc(pZ,dir,type,'bnb',maxOrder)
%    val = supportFunc(pZ,dir,type,'bnbAdv',maxOrder)
%    val = supportFunc(pZ,dir,type,'globOpt',maxOrder,tol)
%    val = supportFunc(pZ,dir,type,'split',splits)
%
% Inputs:
%    pZ - polyZonotope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - minimum ('lower'), maximum ('upper') or range ('range')
%    method - method that is used to calculate the bounds for the dependent
%             part of the polynomial zonotope
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'bnb': taylor models with "branch and bound" algorithm
%              'bnbAdv': taylor models with advandced bnb-algorithm
%              'globOpt': verified global optimization 
%              'bernstein': conversion to a bernstein polynomial
%    maxOrder - maximum polynomial order of the taylor model
%    tol - tolerance for the verified global optimization method
%    split - number of splits that are performed to calculate the bounds
%
% Outputs:
%    val - interval object specifying the upper and lower bound along the
%          direction
%
% Examples:
%    pZ = polyZonotope([0;0],[2 -1 2;0 -2 -3],[],[1 0 1;0 1 3]);
%    dir = [1;2];
%    
%    b1 = supportFunc(pZ,dir,'range','interval');
%    ub1 = conHyperplane(dir,supremum(b1));
%    lb1 = conHyperplane(dir,infimum(b1));
%
%    b2 = supportFunc(pZ,dir,'range','globOpt');
%    ub2 = conHyperplane(dir,supremum(b2));
%    lb2 = conHyperplane(dir,infimum(b2));
%
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(ub1,[1,2],'b');
%    plot(lb1,[1,2],'b');
%    plot(ub2,[1,2],'g');
%    plot(lb2,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      29-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % parse input arguments
    type = 'upper';
    method = 'interval';
    
    if nargin > 2 && ~isempty(varargin{1})
       type = varargin{1}; 
    end
    if nargin > 3 && ~isempty(varargin{2})
       method = varargin{2}; 
    end

    % project the polynomial zonotope onto the direction
    pZ_ = dir' * pZ;
    
    % different methods to calculate the over-approximation
    if strcmp(method,'interval')
       
        val = interval(zonotope(pZ_));
        
        if strcmp(type,'lower')
           val = infimum(val); 
        elseif strcmp(type,'upper')
           val = supremum(val);
        end
        
    elseif strcmp(method,'bernstein')
        
        p = size(pZ_.expMat,1);    
        dom = interval(-ones(p,1),ones(p,1));  

        % dependent generators: convert to bernstein polynomial
        B = poly2bernstein(pZ_.G,pZ_.expMat,dom);

        int1 = interval(min(min(B)),max(max(B)));

        % independent generators: enclose zonotope with interval
        int2 = interval(zonotope([pZ_.c,pZ_.Grest]));

        val = int1 + int2;
        
        if strcmp(type,'lower')
           val = infimum(val); 
        elseif strcmp(type,'upper')
           val = supremum(val);
        end
        
    elseif strcmp(method,'split')
        
        % parse input arguments
        splits = 8;
        if nargin >= 5
           splits = varargin{3}; 
        end
        
        % split the polynomial zonotope multiple times to obtain a better 
        % over-approximation of the real shape
        pZsplit{1} = pZ_;

        for i=1:splits
            qZnew = [];
            for j=1:length(pZsplit)
                res = splitLongestGen(pZsplit{j});
                qZnew{end+1} = res{1};
                qZnew{end+1} = res{2};
            end
            pZsplit = qZnew;
        end
        
        % calculate the interval over-approximation
        Min = inf;
        Max = -inf;
        
        for i = 1:length(pZsplit)
           inter = interval(pZsplit{i});
           Min = min(Min,infimum(inter));
           Max = max(Max,supremum(inter));
        end
        
        % extract the desired bound
        if strcmp(type,'lower')
           val = Min; 
        elseif strcmp(type,'upper')
           val = Max;
        else
           val = interval(Min,Max);
        end
        
    else
        
        % over-approximate the independent generators
        valInd = interval(zonotope([pZ_.c,pZ_.Grest]));
        
        % parse input arguments
        maxOrder = 5;
        if nargin >= 5
           maxOrder = varargin{3}; 
        end

        % determine bounds with taylor models and "branch and bound" or
        % advanced "branch and bound"
        if strcmp(method,'bnb') || strcmp(method,'bnbAdv')
            
            % create taylor models
            n = size(pZ_.expMat,1);
            temp = ones(n,1);
            T = taylm(interval(-temp,temp),maxOrder,[],method);
            
            % calculate bounds of the polynomial part (= dependent
            % generators)
            valDep = interval(polyPart(T,pZ_.G,pZ_.expMat));
            
            % add the two parts from the dependent and independent
            % generators
            val = valDep + valInd;
            
            % extract the desired bound
            if strcmp(type,'lower')
               val = infimum(val); 
            elseif strcmp(type,'upper')
               val = supremum(val);
            end
            
        % determine bounds with global verified optimization
        elseif strcmp(method,'globOpt')          
            
            % parse input arguments
            tol = 1e-3;
            if nargin >= 6
               maxOrder = varargin{3}; 
            end
            
            % remove zero entries from the generator matrix
            ind = find(pZ_.G == 0);
            G = pZ_.G;
            expMat = pZ_.expMat;
            G(:,ind) = [];
            expMat(:,ind) = [];
            
            % remove zero rows from the exponent matrix
            expMat(sum(expMat,2) == 0,:) = [];
            
            % function to opimize
            func = @(x) polyPart(x,G,expMat);
            
            % domain for the optimization variables
            n = size(expMat,1);
            temp = ones(n,1);
            dom = interval(-temp,temp);
            
            % calculate the bounds of the depenedent part
            if strcmp(type,'lower')
                minDep = globVerMinimization(func,dom,tol,maxOrder);
                val = minDep + infimum(valInd);
            elseif strcmp(type,'upper')
                maxDep = globVerMinimization(@(x) -func(x),dom,tol,maxOrder);
                val = -maxDep + supremum(valInd);
            else
                valDep = globVerBounds(func,dom,tol,maxOrder);
                val = valDep + valInd;
            end
            
        else
            error('Wrong value for input argument ''type''');
        end
        
    end
end


% Auxlilary functions -----------------------------------------------------

function val = polyPart(x,G,expMat)

    for i = 1:length(G)
       temp = x(1)^expMat(1,i);
       for j = 2:size(expMat,1)
           temp = temp * x(j)^expMat(j,i);
       end
       if i == 1
           val = G(i) * temp;
       else
           val = val + G(i) * temp;
       end
    end
end


%------------- END OF CODE --------------