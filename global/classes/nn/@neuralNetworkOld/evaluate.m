function res = evaluate(obj,val,varargin)
% evaluate - compute the output of a neural network for the given input
%
% Syntax:  
%    res = evaluate(obj,val)
%    res = evaluate(obj,val,type)
%    res = evaluate(obj,val,type,poly)
%    res = evaluate(obj,val,type,poly,nrGen)
%
% Inputs:
%    obj - object of class neuralNetworkOld
%    val - input represented as a vector or set (class: zonotope)
%    type - type of bounds used for the computation ('exact' or 'approx').
%           Only for classes polyZonotope and conZonotope. The default is 
%           'exact'. 
%    poly - polynomial order for the approximation ('lin', 'quad', or 
%           'cub'). Only for classes polyZonotope and taylm. The default 
%           is 'lin'.
%    nrGen - maximum number of generators. Only for classes zonotope and
%            polyZonotope.
%
% Outputs:
%    res - output of the neural network
%
% References:
%    [1] Singh, G., et al. "Fast and Effective Robustness Certification"
%    [2] Tran, H.-D., et al. "Star-Based Reachability Analysis of Deep 
%        Neural Networks", 2019
%    [3] R. Ivanov, et al. "Verisig 2.0: Verification of Neural Network
%        Controllers Using Taylor Model Preconditioning", 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetworkOld

% Author:       Niklas Kochdumper
% Written:      17-September-2021             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    approx = false;
    if nargin > 2 && ~isempty(varargin{1})
       if strcmp(varargin{1},'approx')
          approx = true;
       elseif ~strcmp(varargin{1},'exact')
           throw(CORAerror('CORA:wrongValue','third',"'approx' or 'exact'"));
       end
    end
    poly = 'lin';
    if nargin > 3 && ~isempty(varargin{2})
        poly = varargin{2};
        if ~ismember(poly,{'lin','quad','cub'})
            throw(CORAerror('CORA:wrongValue','fourth',"'lin','quad' or 'cub'"));
        end
    end
    nrGen = [];
    if nargin > 4 && ~isempty(varargin{3})
        nrGen = varargin{3};
    end

    % different algorithms for different set representations
    if isnumeric(val)
        
        % check correctness of input arguments
        if size(val,1) ~= obj.nrOfInputs
            throw(CORAerror('CORA:wrongValue','second',"a vector or set"));
        end
       
        % loop over all layers
        for i = 1:length(obj.W)
           
            input = obj.W{i}*val + obj.b{i};
            
            switch obj.actFun{i}
                case 'ReLU'
                    val = max(0,input);
                case 'sigmoid'
                    val = exp(input)./(1+exp(input));
                case 'tanh'
                    val = tanh(input);
                case 'identity'
                    val = input;
            end
        end
        
        res = val;
        
    elseif isa(val,'zonotope')
        
        % check correctness of input arguments
        if dim(val) ~= obj.nrOfInputs
            throw(CORAerror('CORA:wrongValue','second',"a vector or set"));
        end
        
        Z = val.Z;
        
        % loop over all layers
        for i = 1:length(obj.W)
            
            % simple case for identity activation functions
            if strcmp(obj.actFun{i},'identity')
               Z = obj.W{i}*Z;
               Z(:,1) = Z(:,1) + obj.b{i};
               continue;
            end
           
            % initialization
            Z_ = zeros(size(obj.W{i},1),size(Z,2));
            d = zeros(size(obj.W{i},1),1);
            
            % loop over all neurons in the current layer
            for j = 1:size(obj.W{i},1)
               
                input = obj.W{i}(j,:)*Z;
                input(:,1) = input(:,1) + obj.b{i}(j);
                
                switch obj.actFun{i}
                    case 'ReLU'
                        [Z_(j,:),d(j)] = ReLU_zono(input);
                    case 'sigmoid'
                        [Z_(j,:),d(j)] = sigmoid_zono(input);
                    case 'tanh'
                        [Z_(j,:),d(j)] = tanh_zono(input);
                end
            end
            
            % order reduction
            if ~isempty(nrGen) && size(Z_,2)-1 > nrGen
               temp = reduce(zonotope(Z_),'girard',nrGen/size(Z_,1));
               Z_ = temp.Z;
            end
            temp = diag(d);
            Z = [Z_ temp(:,d > 0)];
        end
        
        res = zonotope(Z); 
        
    elseif isa(val,'polyZonotope')
       
        % check correctness of input arguments
        if dim(val) ~= obj.nrOfInputs
            throw(CORAerror('CORA:wrongValue','second',"a vector or set"));
        end
        
        % properties of polynomial zonotope
        c = val.c; G = val.G; Grest = val.Grest; 
        expMat = val.expMat; id = val.id; id_ = max(id);
        if isempty(Grest)
            Grest = zeros(size(c));
        end
        ind = find(prod(ones(size(val.expMat))-mod(val.expMat,2),1) == 1);
        ind_ = setdiff(1:size(val.expMat,2),ind);
        
        % linear approximation
        if strcmp(poly,'lin')

            % loop over all layers
            for i = 1:length(obj.W)
                
                % simple case for identity activation functions
                if strcmp(obj.actFun{i},'identity')
                   c = obj.W{i}*c + obj.b{i}; G = obj.W{i}*G; 
                   Grest = obj.W{i}*Grest;
                   continue;
                end
                
                % initialization
                c_ = zeros(size(obj.W{i},1),1);
                G_ = zeros(size(obj.W{i},1),size(G,2));
                Grest_ = zeros(size(obj.W{i},1),size(Grest,2));
                d = zeros(size(obj.W{i},1),1);

                % loop over all neurons in the current layer
                for j = 1:size(obj.W{i},1)

                    c_i = obj.W{i}(j,:)*c + obj.b{i}(j);
                    G_i = obj.W{i}(j,:)*G;
                    Grest_i = obj.W{i}(j,:)*Grest;

                    switch obj.actFun{i}
                        case 'ReLU'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                               ReLU_lin(c_i,G_i,Grest_i,expMat,ind, ...
                                                              ind_,approx);
                        case 'sigmoid'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                               sigmoid_lin(c_i,G_i,Grest_i,expMat,ind, ...
                                                              ind_,approx);
                        case 'tanh'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                               tanh_lin(c_i,G_i,Grest_i,expMat,ind, ...
                                                              ind_,approx);
                    end
                end

                % update properties
                c = c_; G = G_; Grest = Grest_;
                
                % order reduction
                [c,G,Grest,expMat,id,d] = reducePolyZono(c,G,Grest, ...
                                                        expMat,id,d,nrGen);
                temp = diag(d);
                Grest = [Grest, temp(:,d > 0)];
                
                % update indices of all-even exponents (for zonotope encl.)
                if size(G,2) < size(G_,2)
                    ind = find(prod(ones(size(expMat)) - ...
                                                    mod(expMat,2),1) == 1);
                    ind_ = setdiff(1:size(expMat,2),ind);                            
                end
            end

            res = polyZonotope(c,G,Grest,expMat,id);
        
        % quadratic approximation
        elseif strcmp(poly,'quad')
            
            % loop over all layers
            for i = 1:length(obj.W)
                
                % simple case for identity activation functions
                if strcmp(obj.actFun{i},'identity')
                   c = obj.W{i}*c + obj.b{i}; G = obj.W{i}*G; 
                   Grest = obj.W{i}*Grest;
                   continue;
                end

                % linear map defined by the weights and the bias
                c = obj.W{i}*c + obj.b{i};
                G = obj.W{i}*G;
                Grest = obj.W{i}*Grest;

                % for image inputs the input dimension is very large, so we
                % perform an order reduction here to reduce comp. time
                if ~isempty(nrGen) && i == 1 && ...
                                         size(G,2) + size(Grest,2) > nrGen
                    [c,G,Grest,expMat,id] = initialOrderReduction(c,...
                                              G,Grest,expMat,id,id_,nrGen);
                    ind = find(prod(ones(size(expMat)) - ...
                                                    mod(expMat,2),1) == 1);
                    ind_ = setdiff(1:size(expMat,2),ind);
                end

                % initialization
                h = size(G,2); q = size(Grest,2);
                temp = q + q*h + 0.5*(q^2 + q);
                c_ = zeros(size(obj.W{i},1),1);
                G_ = zeros(size(obj.W{i},1),0.5*(h^2+h)+h);
                Grest_ = zeros(size(obj.W{i},1),temp);
                d = zeros(size(obj.W{i},1),1);

                % loop over all neurons in the current layer
                for j = 1:size(obj.W{i},1)

                    switch obj.actFun{i}
                        case 'ReLU'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                                ReLU_quad(c(j),G(j,:),Grest(j,:),expMat,...
                                                          ind,ind_,approx);
                        case 'sigmoid'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                                sigmoid_quad(c(j),G(j,:),Grest(j,:), ...
                                                   expMat,ind,ind_,approx);
                        case 'tanh'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                              tanh_quad(c(j),G(j,:),Grest(j,:),expMat, ...
                                                          ind,ind_,approx);
                    end
                end

                % update properties
                c = c_; G = G_; Grest = Grest_(:,sum(abs(Grest_),1) > 0);
                expMat = [expMat,squareExpMat(expMat)];
                [expMat,G] = removeRedundantExponents(expMat,G);
                
                % order reduction
                [c,G,Grest,expMat,id,d] = reducePolyZono(c,G,Grest, ...
                                                      expMat,id,d,nrGen);
                temp = diag(d);
                Grest = [Grest, temp(:,d > 0)];
                
                % update indices of all-even exponents (for zonotope encl.)
                ind = find(prod(ones(size(expMat))-mod(expMat,2),1) == 1);
                ind_ = setdiff(1:size(expMat,2),ind);
            end

            res = polyZonotope(c,G,Grest,expMat,id);
            
        % cubic approximation
        elseif strcmp(poly,'cub')
            
            % loop over all layers
            for i = 1:length(obj.W)
                
                % simple case for identity activation functions
                if strcmp(obj.actFun{i},'identity')
                   c = obj.W{i}*c + obj.b{i}; G = obj.W{i}*G; 
                   Grest = obj.W{i}*Grest;
                   continue;
                end

                % linear map defined by the weights and the bias
                c = obj.W{i}*c + obj.b{i};
                G = obj.W{i}*G;
                Grest = obj.W{i}*Grest;

                % for image inputs the input dimension is very large, so we
                % perform an order reduction here to reduce comp. time
                if ~isempty(nrGen) && i == 1 && ...
                                          size(G,2) + size(Grest,2) > nrGen
                    [c,G,Grest,expMat,id] = initialOrderReduction(c,...
                                              G,Grest,expMat,id,id_,nrGen);
                    ind = find(prod(ones(size(expMat)) - ...
                                                    mod(expMat,2),1) == 1);
                    ind_ = setdiff(1:size(expMat,2),ind);
                end

                % initialization
                h = size(G,2); q = size(Grest,2);
                temp = nchoosek(3+q,q)-1 + h*q + ...
                                        0.5*(h^2+h)*q + 0.5*(q^2+q)*h;
                c_ = zeros(size(obj.W{i},1),1);
                G_ = zeros(size(obj.W{i},1),nchoosek(3+h,h)-1);
                Grest_ = zeros(size(obj.W{i},1),temp);
                d = zeros(size(obj.W{i},1),1);

                % loop over all neurons in the current layer
                for j = 1:size(obj.W{i},1)

                    switch obj.actFun{i}
                        case 'ReLU'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                               ReLU_cub(c(j),G(j,:),Grest(j,:),expMat, ...
                                                          ind,ind_,approx);
                        case 'sigmoid'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                            sigmoid_cub(c(j),G(j,:),Grest(j,:),expMat, ...
                                                          ind,ind_,approx);
                        case 'tanh'
                            [c_(j),G_(j,:),Grest_(j,:),d(j)] = ...
                               tanh_cub(c(j),G(j,:),Grest(j,:),expMat, ...
                                                          ind,ind_,approx);
                    end
                end

                % update properties
                c = c_; G = G_; Grest = Grest_(:,sum(abs(Grest_),1) > 0);  
                expMat = [expMat,squareExpMat(expMat),cubeExpMat(expMat)];
                [expMat,G] = removeRedundantExponents(expMat,G);
                
                % order reduction
                [c,G,Grest,expMat,id,d] = reducePolyZono(c,G,Grest, ...
                                                      expMat,id,d,nrGen);
                temp = diag(d);
                Grest = [Grest, temp(:,d > 0)];
                
                % update indices of all-even exponents (for zonotope encl.)
                ind = find(prod(ones(size(expMat))-mod(expMat,2),1) == 1);
                ind_ = setdiff(1:size(expMat,2),ind);
            end

            res = polyZonotope(c,G,Grest,expMat,id);
        end
        
    elseif isa(val,'taylm')
        
        % loop over all layers
        for i = 1:length(obj.W)
            
            % linear transformation defined by the weights and biases
            val = obj.W{i}*val + obj.b{i};
            
            % simple case for identity activation functions
            if strcmp(obj.actFun{i},'identity')
               continue;
            end
            
            % loop over all neurons in the current layer
            for j = 1:size(val,1)
                
                switch obj.actFun{i}
                    case 'ReLU'
                        val(j) = ReLU_taylm(val(j),poly);
                    case 'sigmoid'
                        val(j) = sigmoid_taylm(val(j),poly);
                    case 'tanh'
                        val(j) = tanh_taylm(val(j),poly); 
                end
            end
        end
        
        res = val;
        
    elseif isa(val,'conZonotope')
        
        % check activation function
        temp = cellfun(@(x) ismember(x,{'ReLU','identity'}),obj.actFun);
        if ~all(temp)
           ind = find(temp == 0);
           throw(CORAerror('CORA:notSupported',...
               ['Activation function "',obj.actFun{ind(1)}, ...
                  '" is not supported for conZonotope objects!'])); 
        end
        
        % convert constrained zonotope to star set
        [c,G,C,d,l,u] = conversionConZonoStarSet(val);
        
        % predefine options for linear programming for speed-up
        options = optimoptions('linprog','display','off');
                
        % loop over all layers
        for i = 1:length(obj.W)
            
            % linear transformation defined by the weights and biases
            c = obj.W{i}*c + obj.b{i}; G = obj.W{i}*G;
            
            % simple case for identity activation functions
            if strcmp(obj.actFun{i},'identity')
               continue;
            end
            
            % loop over all neurons in the current layer
            for j = 1:length(c)
                [c,G,C,d,l,u] = ReLU_starSet(c,G,C,d,l,u,j,options,approx);
            end            
        end
        
        % convert star set to constrained zonotope
        res = conversionStarSetConZono(c,G,C,d,l,u);
            
    else
        throw(CORAerror('CORA:notSupported','Given set representation not supported.'));
    end
end


% Auxiliary Functions -----------------------------------------------------

function [res,d] = ReLU_zono(input)
% enclose the ReLU activation function with a zonotope according to 
% Theorem 3.1 in [1]

    % compute lower and upper bound
    l = input(1) - sum(abs(input(2:end)));
    u = input(1) + sum(abs(input(2:end)));

    % compute output
    if u < 0
        res = zeros(size(input)); d = 0;
    elseif l > 0
        res = input; d = 0;
    else
        lambda = u/(u-l); mu = -(u*l)/(2*(u-l));
        res = lambda*input;
        res(1) = res(1) + mu;
        d = mu;
    end
end

function [res,d] = sigmoid_zono(input)
% enclose the sigmoid activation function with a zonotope according to 
% Theorem 3.2 in [1]

    % sigmoid function and its derivative
    f = @(x) exp(x)/(1+exp(x));
    df = @(x) f(x)*(1-f(x));
    
    % compute lower and upper bound
    l = input(1) - sum(abs(input(2:end)));
    u = input(1) + sum(abs(input(2:end)));

    % compute auxiliary functions
    lambda = min(df(l),df(u));
    mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
    mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
    
    % compute output
    if l == u
        res = zeros(1,length(input));
        res(1) = f(u); d = 0;
    else
        res = lambda*input;
        res(1) = res(1) + mu1;
        d = mu2;
    end
end

function [res,d] = tanh_zono(input)
% enclose the hyperbolic tangent activation function with a zonotope 
% according to Theorem 3.2 in [1]

    % hyperbolic tangent function and its derivative
    f = @(x) tanh(x);
    df = @(x) 1 - tanh(x)^2;
    
    % compute lower and upper bound
    l = input(1) - sum(abs(input(2:end)));
    u = input(1) + sum(abs(input(2:end)));

    % compute auxiliary functions
    lambda = min(df(l),df(u));
    mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
    mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
    
    % compute output
    if l == u
        res = zeros(1,length(input));
        res(1) = f(u); d = 0;
    else
        res = lambda*input;
        res(1) = res(1) + mu1;
        d = mu2;
    end
end

function [c,G,Grest,d] = ReLU_lin(c,G,Grest,E,ind,ind_,approx)
% enclose the ReLU activation function with a polynomial zonotope according 
% to Theorem 3.1 in [1]

    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);

    % compute output
    if u < 0
        c = 0; G = zeros(size(G)); 
        Grest = zeros(size(Grest)); d = 0;
    elseif l > 0
        d = 0;
    else
        lambda = u/(u-l); mu = -(u*l)/(2*(u-l));
        c = lambda*c + mu; G = lambda*G;
        Grest = lambda*Grest;
        d = mu;
    end
end

function [c,G,Grest,d] = sigmoid_lin(c,G,Grest,E,ind,ind_,approx)
% enclose the sigmoid activation function with a polynomial zonotope 
% according to Theorem 3.2 in [1]

    % sigmoid function and its derivative
    f = @(x) exp(x)/(1+exp(x));
    df = @(x) f(x)*(1-f(x));

    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);

    % compute auxiliary functions
    lambda = min(df(l),df(u));
    mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
    mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
    
    % compute output
    if l == u
        G = zeros(size(G)); Grest = zeros(size(Grest));
        c = f(u); d = 0;
    else
        c = lambda*c + mu1; G = lambda*G;
        Grest = lambda*Grest;
        d = mu2;
    end
end

function [c,G,Grest,d] = tanh_lin(c,G,Grest,E,ind,ind_,approx)
% enclose the hyperbolic tangent activation function with a polynomial 
% zonotope according to Theorem 3.2 in [1]

    % hyperbolic tangent function and its derivative
    f = @(x) tanh(x);
    df = @(x) 1 - tanh(x)^2;

    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);

    % compute auxiliary functions
    lambda = min(df(l),df(u));
    mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
    mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
    
    % compute output
    if l == u
        G = zeros(size(G)); Grest = zeros(size(Grest));
        c = f(u); d = 0;
    else
        c = lambda*c + mu1; G = lambda*G;
        Grest = lambda*Grest;
        d = mu2;
    end
end

function [c,G,Grest,d] = ReLU_quad(c,G,Grest,E,ind,ind_,approx)
% enclose the ReLU activation function with a polynomial zonotope by
% fitting a quadratic function

    % properties
    h = length(G); q = length(Grest);

    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute output
    if u < 0
        c = 0; G = zeros(1,h + 0.5*(h^2 + h)); 
        Grest = zeros(1,q + q*h + 0.5*(q^2 + q));
        d = 0;
    elseif l > 0
        G = [G, zeros(1,0.5*(h^2 + h))];
        Grest = [Grest,zeros(1,q*h + 0.5*(q^2 + q))];
        d = 0;
    else
        
        % compute quadratic function that best fits the activation function
        c_a = u/(u-l)^2; c_b = -2*l*u/(u-l)^2; 
        c_c = u^2*(2*l-u)/(u-l)^2 + u;

        % compute difference between ReLU and quadratic approximation
        L1 = minMaxQuadFun(-c_a,1-c_b,-c_c,0,u);
        L2 = minMaxQuadFun(-c_a,-c_b,-c_c,l,0);
        L = L1 | L2; d = rad(L);

        % evaluate the quadratic approximation on the polynomial zonotope
        [c,G,Grest] = quadApproxPolyZono(c,G,Grest,c_a,c_b,c_c);
        c = c + center(L);
    end
end

function [c,G,Grest,d] = sigmoid_quad(c,G,Grest,E,ind,ind_,approx)
% enclose the sigmoid activation function with a polynomial zonotope by
% fitting a quadratic function

    % sigmoid activation function
    f = @(x) exp(x)./(1+exp(x));
    
    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute quadratic function that best fits the activation function
    x = linspace(l,u,10); y = f(x);
    
    temp = leastSquarePolyFunc(x,y,2);
    c_c = temp(1); c_b = temp(2); c_a = temp(3);
    
    % compute the difference between activation function and quad. fit    
    L = minMaxDiffQuad(c_a,c_b,c_c,l,u,'sigmoid');
    d = rad(L);

    % evaluate the quadratic approximation on the polynomial zonotope
    [c,G,Grest] = quadApproxPolyZono(c,G,Grest,c_a,c_b,c_c);
    c = c + center(L);
end

function [c,G,Grest,d] = tanh_quad(c,G,Grest,E,ind,ind_,approx)
% enclose the hyperbolic tangent activation function with a polynomial 
% zonotope by fitting a quadratic function

    % hyperbolic tangent activation function
    f = @(x) tanh(x);
    
    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute quadratic function that best fits the activation function
    x = linspace(l,u,10); y = f(x);
    temp = leastSquarePolyFunc(x,y,2);
    c_c = temp(1); c_b = temp(2); c_a = temp(3);
    
    % compute the difference between activation function and quad. fit
    L = minMaxDiffQuad(c_a,c_b,c_c,l,u,'tanh');
    d = rad(L);
    
    % evaluate the quadratic approximation on the polynomial zonotope
    [c,G,Grest] = quadApproxPolyZono(c,G,Grest,c_a,c_b,c_c);     
    c = c + center(L);
end

function [c,G,Grest,d] = ReLU_cub(c,G,Grest,E,ind,ind_,approx)
% enclose the ReLU activation function with a polynomial zonotope by
% fitting a cubic function

    % properties
    h = length(G); q = length(Grest);

    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute output
    if u < 0
        c = 0; G = zeros(1,nchoosek(3+h,h)-1); 
        temp = nchoosek(3+q,q)-1 + h*q + 0.5*(h^2+h)*q + 0.5*(q^2+q)*h;
        Grest = zeros(1,temp); d = 0;
    elseif l > 0
        G = [G, zeros(1,nchoosek(3+h,h)-1-h)];
        temp = nchoosek(3+q,q)-1 + h*q + 0.5*(h^2+h)*q + 0.5*(q^2+q)*h-q;
        Grest = [Grest,zeros(1,temp)]; d = 0; 
    else
        
        % compute quadratic function that best fits the activation function
        x = linspace(l,u,10); y = max(0*x,x);
        temp =  leastSquarePolyFunc(x,y,3);
        c_d = temp(1); c_c = temp(2); c_b = temp(3); c_a = temp(4);

        % compute difference between ReLU and cubic approximation
        L1 = minMaxCubFun(-c_a,-c_b,1-c_c,-c_d,0,u);
        L2 = minMaxCubFun(-c_a,-c_b,-c_c,-c_d,l,0);
        L = L1 | L2; d = rad(L);

        % evaluate the cubic approximation on the polynomial zonotope
        [c,G,Grest] = cubApproxPolyZono(c,G,Grest,c_a,c_b,c_c,c_d);
        c = c + center(L);
    end
end

function [c,G,Grest,d] = sigmoid_cub(c,G,Grest,E,ind,ind_,approx)
% enclose the sigmoid activation function with a polynomial zonotope by
% fitting a cubic function

    % sigmoid activation function
    f = @(x) exp(x)./(1+exp(x));
    
    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute cubic function that best fits the activation function
    x = linspace(l,u,10); y = f(x);
    temp =  leastSquarePolyFunc(x,y,3);
    c_d = temp(1); c_c = temp(2); c_b = temp(3); c_a = temp(4);
    
    % compute the difference between activation function and cubic fit    
    L = minMaxDiffCub(c_a,c_b,c_c,c_d,l,u,'sigmoid');
    d = rad(L);
    
    % evaluate the cubic approximation on the polynomial zonotope
    [c,G,Grest] = cubApproxPolyZono(c,G,Grest,c_a,c_b,c_c,c_d);
    c = c + center(L);
end

function [c,G,Grest,d] = tanh_cub(c,G,Grest,E,ind,ind_,approx)
% enclose the hyperbolic tangent activation function with a polynomial 
% zonotope by fitting a cubic function

    % sigmoid activation function
    f = @(x) tanh(x);
    
    % compute lower and upper bound
    [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx);
    
    % compute cubic function that best fits the activation function
    x = linspace(l,u,10); y = f(x);
    temp =  leastSquarePolyFunc(x,y,3);
    c_d = temp(1); c_c = temp(2); c_b = temp(3); c_a = temp(4);
    
    % compute the difference between activation function and cubic fit    
    L = minMaxDiffCub(c_a,c_b,c_c,c_d,l,u,'tanh');
    d = rad(L);

    % evaluate the cubic approximation on the polynomial zonotope
    [c,G,Grest] = cubApproxPolyZono(c,G,Grest,c_a,c_b,c_c,c_d);
    c = c + center(L);
end

function res = ReLU_taylm(tay,poly)
% enclose the ReLU activation function with a Taylor model by
% fitting a quadratic function

    % compute lower and upper bound
    int = interval(tay);
    l = infimum(int); u = supremum(int);
    
    % compute output
    if u < 0
        res = 0*tay;
    elseif l > 0
        res = tay; 
    else
        
        if strcmp(poly,'lin')
            
            % compute linear enclosure according to Theorem 3.1 in [1]
            lambda = u/(u-l); mu = -(u*l)/(2*(u-l));
            res = lambda*tay + interval(0,2*mu);           
            
        elseif strcmp(poly,'quad')
            
            % compute quadratic fit to the activation function
            x = linspace(l,u,10); y = max(0*x,x);
            temp = leastSquarePolyFunc(x,y,2);
            c_c = temp(1); c_b = temp(2); c_a = temp(3);

            % compute difference between ReLU and quadratic approximation
            L1 = minMaxQuadFun(-c_a,1-c_b,-c_c,0,u);
            L2 = minMaxQuadFun(-c_a,-c_b,-c_c,l,0);
            L = L1 | L2;

            % compute resulting Taylor model
            res = c_a*tay^2 + c_b*tay + c_c + L;
            
        elseif strcmp(poly,'cub')
            
            % compute cubic fit to the activation function
            x = linspace(l,u,10); y = max(0*x,x);
            temp =  leastSquarePolyFunc(x,y,3);
            c_d = temp(1); c_c = temp(2); c_b = temp(3); c_a = temp(4);

            % compute difference between ReLU and cubic approximation
            L1 = minMaxCubFun(-c_a,-c_b,1-c_c,-c_d,0,u);
            L2 = minMaxCubFun(-c_a,-c_b,-c_c,-c_d,l,0);
            L = L1 | L2;

            % compute resulting Taylor model
            res = c_a*tay^3 + c_b*tay^2 + c_c*tay + c_d + L;
        end
    end
end

function res = sigmoid_taylm(tay,poly)
% enclose the sigmoid activation function with a Taylor model

    % sigmoid activation function and its derivatives
    f = @(x) exp(x)./(1+exp(x));
    df = @(x) f(x)*(1-f(x));
    df2 = @(x) f(x)*(1-f(x))*(1-2*f(x));
    df3 = @(x) f(x)*(1-f(x))*(1-6*f(x)*(1-f(x)));
    
    % compute lower and upper bound
    int = interval(tay);
    l = infimum(int); u = supremum(int); c = center(int);
    
    % consider different polynomial orders
    if strcmp(poly,'lin')
        
        % compute linear enclosure according to Theorem 3.2 in [1]
        lambda = min(df(l),df(u));
        mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
        mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
        
        res = lambda*tay + interval(mu1-mu2,mu1+mu2);
        
    elseif strcmp(poly,'quad')
    
        % obtain quadratic approximation using the Taylor series expansion
        % according to Example 1 in [3]
        c_a = 0.5*df2(c);
        c_b = df(c) - df2(c)*c;
        c_c = f(c) - df(c)*c + 0.5*df2(c)*c^2;

        % compute the difference between activation function and quad. fit    
        L = minMaxDiffQuad(c_a,c_b,c_c,l,u,'sigmoid');

        % compute resulting Taylor model
        res = c_a * tay^2 + c_b * tay + c_c + L;
    
    elseif strcmp(poly,'cub')
        
        % compute cubic approximation using the Taylor series expansion
        c_a = 1/6*df3(c);
        c_b = 0.5*df2(c) - 0.5*df3(c)*c;
        c_c = df(c) - df2(c)*c + 0.5*c^2*df3(c);
        c_d = f(c) - df(c)*c + 0.5*df2(c)*c^2 - 1/6*df3(c)*c^3;

        % compute the difference between activation function and cubic fit    
        L = minMaxDiffCub(c_a,c_b,c_c,c_d,l,u,'sigmoid');

        % compute resulting Taylor model
        res = c_a*tay^3 + c_b*tay^2 + c_c*tay + c_d + L;
    end
end

function res = tanh_taylm(tay,poly)
% enclose the hyperbolic tangent activation function with a Taylor model

    % sigmoid activation function
    f = @(x) tanh(x);
    df = @(x) 1 - tanh(x)^2;
    df2 = @(x) -2*tanh(x)*(1-tanh(x)^2); 
    df3 = @(x) -2*(tanh(x)^2-1)^2 - 4*tanh(x)^2*(tanh(x)^2-1);
    
    % compute lower and upper bound
    int = interval(tay);
    l = infimum(int); u = supremum(int); c = center(int);
    
    % consider different polynomial orders
    if strcmp(poly,'lin')
    
        % compute linear enclosure according to Theorem 3.2 in [1]
        lambda = min(df(l),df(u));
        mu1 = 0.5*(f(u) + f(l) - lambda*(u+l));
        mu2 = 0.5*(f(u) - f(l) - lambda*(u-l));
        
        res = lambda*tay + interval(mu1-mu2,mu1+mu2);
    
    elseif strcmp(poly,'quad')
    
        % obtain quadratic approximation using the Taylor series expansion
        % according to Example 1 in [3]
        c_a = 0.5*df2(c);
        c_b = df(c) - df2(c)*c;
        c_c = f(c) - df(c)*c + 0.5*df2(c)*c^2;

        % compute the difference between activation function and quad. fit    
        L = minMaxDiffQuad(c_a,c_b,c_c,l,u,'tanh');

        % compute resulting Taylor model
        res = c_a * tay^2 + c_b * tay + c_c + L;
    
    elseif strcmp(poly,'cub')
        
        % compute cubic approximation using the Taylor series expansion
        c_a = 1/6*df3(c);
        c_b = 0.5*df2(c) - 0.5*df3(c)*c;
        c_c = df(c) - df2(c)*c + 0.5*c^2*df3(c);
        c_d = f(c) - df(c)*c + 0.5*df2(c)*c^2 - 1/6*df3(c)*c^3;

        % compute the difference between activation function and cubic fit    
        L = minMaxDiffCub(c_a,c_b,c_c,c_d,l,u,'tanh');

        % compute resulting Taylor model
        res = c_a*tay^3 + c_b*tay^2 + c_c*tay + c_d + L;
    end
end

function res = ReLU_conZono(cZ,i)
% enclose the ReLU activation function with a constrained zonotope based on
% the results for star sets in [2]
    
    n = dim(cZ);

    % get lower and upper bound
    int = interval(project(cZ,i));
    l = infimum(int); u = supremum(int);

    % compute output according to Sec. 3.2 in [2]
    if u < 0
        M = eye(n); M(:,i) = zeros(n,1);
        res = M * cZ;
    elseif l > 0
        res = cZ;
    else        
        c = cZ.Z(:,1); G = cZ.Z(:,2:end);
        a = u*l/(u-l);
        
        A = [-G(i,:) u/2 -l/2 0; (-u/(u-l))*G(i,:) u/2 0 a/2];
        
        if ~isempty(cZ.A)
            A = [cZ.A zeros(size(cZ.A,1),3); A];
        end
        
        b = [cZ.b; -(u+l)/2+c(i); -u/2-a/2+u*c(i)/(u-l)];
        c(i) = u/2; 
        G = [G, zeros(size(G,1),3)];
        G(i,:) = [zeros(1,size(G,2)-3), u/2, 0, 0];
        
        res = conZonotope(c,G,A,b);
    end
end

function [c,G,C,d,l_,u_] = ReLU_starSet(c,G,C,d,l_,u_,i,options,approx)
% enclose the ReLU activation function with a constrained zonotope based on
% the results for star sets in [2]
    
    n = length(c); m = size(G,2); M = eye(n); M(:,i) = zeros(n,1);

    % get lower bound
    if approx
        c_ = c(i) + 0.5*G(i,:)*(u_-l_); G_ = 0.5*G(i,:)*diag(u_-l_); 
        l = c_ - sum(abs(G_));
    else
        [~,temp] = linprog(G(i,:),C,d,[],[],[],[],options);
        l = c(i) + temp;
    end

    % compute output according to Sec. 3.2 in [2]
    if l < 0

        % compute upper bound
        if approx
            u = c_ + sum(abs(G_));
        else
            [~,temp] = linprog(-G(i,:),C,d,[],[],[],[],options);
            u = c(i) - temp;
        end
        
        if u <= 0
            c = M * c; G = M * G;
        else
            C1 = [zeros(1,m),-1]; d1 = 0;
            C2 = [G(i,:) -1]; d2 = -c(i);
            C3 = [-u/(u-l)*G(i,:) 1]; d3 = u/(u-l)*(c(i)-l);
            C0 = [C zeros(size(C,1),1)]; d0 = d;
            C = [C0; C1; C2; C3]; d = [d0; d1; d2; d3];
            temp = zeros(n,1); temp(i) = 1;
            c = M*c; G = M*G; G = [G, temp];
            l_ = [l_;0]; u_ = [u_;u];
        end
    end
end

function coeff = leastSquarePolyFunc(x,y,n)
% determine the optimal cubic function fit that minimizes the squared
% distance to the data points

    A = x'.^(0:n);
    coeff = pinv(A)*y';
end

function G_ = squareGenMat(G)
% compute the generator matrix for the squared polynomial zonotope

    if ~isempty(G)
        n = length(G);
        G_ = zeros(1,0.5*n*(n+1));
        temp = G'*G; cnt = n;

        for i = 1:n-1
           G_(i) = temp(i,i);
           G_(cnt+1:cnt+n-i) = 2*temp(i,i+1:n);
           cnt = cnt + n-i;
        end
        G_(n) = temp(end,end);
    else
       G_ = []; 
    end
end

function E = squareExpMat(E)
% compute the exponent matrix for the squared polynomial zonotope

    E_ = [];
    
    for i = 1:size(E,2)-1
        for j = i+1:size(E,2)
           E_ = [E_ E(:,i)+E(:,j)];
        end
    end

    E = [2*E E_];
end

function G = cubeGenMat(G)
% compute the generator matrix for the cubed polynomial zonotope

    G1 = []; G2 = []; G3 = [];
    
    for i = 1:size(G,2)-1
       for j = i+1:size(G,2)
          G1 = [G1 3*G(i)^2*G(j)]; 
       end
    end
    
    for i = 1:size(G,2)
       for j = 1:i-1
          G2 = [G2 3*G(i)^2*G(j)]; 
       end
    end
    
    for i = 1:size(G,2)-1
       for j = i+1:size(G,2)-1
           for k = j+1:size(G,2)
               G3 = [G3 6*G(i)*G(j)*G(k)];
           end
       end
    end
    
    G = [G.^3 G1 G2 G3];
end

function E = cubeExpMat(E)
% compute the exponent matrix for the cubed polynomial zonotope

    E1 = []; E2 = []; E3 = [];
    
    for i = 1:size(E,2)-1
       for j = i+1:size(E,2)
          E1 = [E1 2*E(:,i)+E(:,j)]; 
       end
    end
    
    for i = 2:size(E,2)
       for j = 1:i-1
          E2 = [E2 2*E(:,i)+E(:,j)]; 
       end
    end
    
    for i = 1:size(E,2)-1
       for j = i+1:size(E,2)-1
           for k = j+1:size(E,2)
               E3 = [E3 E(:,i)+E(:,j)+E(:,k)];
           end
       end
    end
    
    E = [3*E E1 E2 E3];
end

function int = minMaxQuadFun(a,b,c,l,u)
% compute the maximimum and minimum value of a quadratic function on an
% interval

   % quadratic function
   f = @(x) a*x^2 + b*x + c;

   % extremal point of the function
   x_ = -b/(2*a);
   
   % maximum and minimum
   if x_ > l && x_ < u
      val = [f(l),f(x_),f(u)];
   else
      val = [f(l),f(u)];
   end

   int = interval(min(val),max(val));
end

function int = minMaxCubFun(a,b,c,d,l,u)
% compute the maximimum and minimum value of a cubic function on an
% interval

   % cubic function
   f = @(x) a*x^3 + b*x^2 + c*x + d;

   % potential extremal points of the function
   x_ = roots([3*a 2*b c]);
   
   % maximum and minimum
   val = [f(l),f(u)];
   for i = 1:length(x_)
       if isreal(x_(i)) && x_(i) > l && x_(i) < u
          val = [val,f(x_(i))];
       end
   end
   int = interval(min(val),max(val));
end

function int = minMaxDiffQuad(a,b,c,l,u,type)
% compute the maximum and the minimum difference between the activation
% function and a quadratic fit

    tol = 0.001;
    
    % compute rough bounds for the derivative
    if strcmp(type,'tanh')
        der1 = interval(0,1);                       % hyperbolic tangent
    elseif strcmp(type,'sigmoid')
        der1 = interval(0,0.25);                    % sigmoid
    end
    
    der2 = -(b + 2*a*interval(l,u));                % quadratic fit
    
    der = supremum(abs(der1 - der2));

    % determine function bounds by sampling
    dx = tol/der;
    x = linspace(l,u,ceil((u-l)/dx));
    
    if strcmp(type,'tanh')
        y = tanh(x) - a*x.^2 - b*x - c;
    elseif strcmp(type,'sigmoid')
        y = 1./(1+exp(-x)) - a*x.^2 - b*x - c;
    end
    
    int = interval(min(y)-tol,max(y)+tol);
end

function int = minMaxDiffCub(a,b,c,d,l,u,type)
% compute the maximum and the minimum difference between the activation
% function and a cubic fit

    tol = 0.001;
    
    % compute rough bounds for the derivative
    if strcmp(type,'tanh')
        der1 = interval(0,1);                       % hyperbolic tangent
    elseif strcmp(type,'sigmoid')
        der1 = interval(0,0.25);                    % sigmoid
    end
    
    der2 = minMaxQuadFun(3*a,2*b,c,l,u);            % cubic fit
    
    der = supremum(abs(der1 - der2));

    % determine function bounds by sampling
    dx = tol/der;
    x = linspace(l,u,ceil((u-l)/dx));
    
    if strcmp(type,'tanh')
        y = tanh(x) - a*x.^3 - b*x.^2 - c*x - d;
    elseif strcmp(type,'sigmoid')
        y = 1./(1+exp(-x)) - a*x.^3 - b*x.^2 - c*x - d;
    end
    
    int = interval(min(y)-tol,max(y)+tol);
end

function [c,G,C,d,l,u] = conversionConZonoStarSet(cZ)
% convert a constrained zonotope to a star set
    
    m = size(cZ.Z,2)-1;

    if isempty(cZ.A)
        
        c = cZ.Z(:,1); G = cZ.Z(:,2:end);
        C = [eye(m); -eye(m)]; d = ones(2*m,1);
        u = ones(m,1); l = -ones(m,1);
        
    else
        % compute point satisfying all constraints with pseudo inverse
        p_ = pinv(cZ.A)*cZ.b;

        % compute null-space of constraints
        T = null(cZ.A);

        % transform boundary constraints of the factor hypercube
        A = [eye(m);-eye(m)]; b = ones(2*m,1);
        C = A*T; d = b - A*p_;
        c = cZ.Z(:,1) + cZ.Z(:,2:end)*p_; G = cZ.Z(:,2:end)*T;
        int = interval(mptPolytope(C,d));
        l = infimum(int); u = supremum(int);
    end
end

function cZ = conversionStarSetConZono(c,G,C,d,l_,u_)
% convert a star set to a constrained zonotope  

    % normalize halfspace normal vector length
    temp = sqrt(sum(C.^2,2)); C = diag(temp)*C; d = temp.*d;

    m = size(G,2);

    % find hypercube constraints
    u = zeros(m,1); ind = find(any(C' == 1)); 
    ind_ = setdiff(1:size(C,1),ind); C_ = C(ind,:); d_ = d(ind);
    
    for i = 1:m
       ind = find(C_(:,i) == 1);
       if isempty(ind)
          u(i) = u_(i); 
       else
          u(i) = min(d_(ind)); 
       end
    end
    
    C = C(ind_,:); d = d(ind_);

    l = zeros(m,1); ind = find(any(C' == -1)); 
    ind_ = setdiff(1:size(C,1),ind); C_ = C(ind,:); d_ = d(ind);
    
    for i = 1:m
       ind = find(C_(:,i) == -1);
       if isempty(ind)
          l(i) = l_(i);
       else
          l(i) = max(-d_(ind)); 
       end
    end
    
    C = C(ind_,:); d = d(ind_);
    
    % scale constraints according to hypercube dimensions
    int = interval(l,u); cen = center(int); R = diag(rad(int));
    d = d - C*cen; C = C*R; c = c + G*cen; G = G*R;
    
    % represent inequality constraints as equivalent equality constraints
    if ~isempty(C)
        Z = zonotope(interval(-ones(m,1),ones(m,1)));
        l = infimum(interval(C*Z));
        int = interval(l,d); cen = center(int); r = rad(int);
        A = [C,diag(r)]; b = cen; G = [G, zeros(size(G,1),length(r))];
    else
       A = []; b = []; 
    end
    
    % construct resulting constrained zonotope
    cZ = conZonotope(c,G,A,b);
end

function [l,u] = compBoundsPolyZono(c,G,Grest,E,ind,ind_,approx)
% compute the lower and upper bound of a polynomial zonotope

    if approx
        c_ = c + 0.5*sum(G(:,ind),2);
    
        l = c_-sum(abs(0.5*G(:,ind)))-sum(abs(G(:,ind_)))-sum(abs(Grest));
        u = c_+sum(abs(0.5*G(:,ind)))+sum(abs(G(:,ind_)))+sum(abs(Grest));
    else
        pZ = polyZonotope(c,G,Grest,E);
        int = interval(pZ,'split');
        l = infimum(int); u = supremum(int);
    end
end

function [c,G,Grest] = quadApproxPolyZono(c,G,Grest,c_a,c_b,c_c)
% evaluate a quadratic function c_a*x^2 + c_b*x + c_c on a polynomial
% zonotope
        
    q = length(Grest);

    % compute (c + G + Grest)*(c + G + Grest) = 
    %         = c*c + 2*c*G + G*G + 2*c*Grest + 2*G*Grest + Grest*Grest
    G_quad = squareGenMat(G);
    Grest_quad = squareGenMat(Grest);
    Grest_quad(:,1:q) = 0.5*Grest_quad(:,1:q);

    % construct resulting polynomial zonotope
    Grest = [c_b*Grest + 2*c_a*c*Grest, c_a*Grest_quad, ...
                                            2*c_a*reshape(G'*Grest,1,[])];
    G = [c_b*G + 2*c_a*c*G, c_a*G_quad];
    c = c_a*(c^2 + sum(Grest_quad(:,1:q))) + c_b*c + c_c;
end

function [c,G,Grest] = cubApproxPolyZono(c,G,Grest,c_a,c_b,c_c,c_d)
% evaluate a cubic function c_a*x^3 + c_b*x^2 + c_c*x + c_d on a polynomial
% zonotope

    h = length(G); q = length(Grest);
    
    % compute (c + G + Grest)^3 = (c^3 + 3*c^2*G + 3*c*G*G + G*G*G) + ...
    %  3*c^2*Grest + 6*c*G*Grest + 3*c*Grest*Grest + 3*G*Grest*G + ...
    %  3*Grest*Grest*G + G_r*G_r*G_r
    G_quad = squareGenMat(G);
    G_quad_ = G_quad; G_quad_(:,1:h) = 0.5*G_quad_(:,1:h);
    c_quad = 0.5*sum(G_quad(:,1:h));
    G_cub = cubeGenMat(G);
    Grest_quad = squareGenMat(Grest);
    Grest_quad(:,1:q) = 0.5*Grest_quad(:,1:q);
    crest_quad = 0.5*sum(Grest_quad(:,1:q));
    Grest_cub = cubeGenMat(Grest);

    % construct resulting polynomial zonotope
    Grest = [(c_a*3*c^2 + c_b*2*c + c_c + c_a*3*c_quad)*Grest, ...
             (c_a*6*c + c_b*2)*reshape(G'*Grest,1,[]), ...
             (c_a*3*c + c_b)*Grest_quad, ...
             (c_a*3)*reshape(Grest'*G_quad_,1,[]), ...
             (c_a*3)*reshape(G'*Grest_quad,1,[]), c_a*Grest_cub];
    G = [(c_a*3*c^2 + c_b*2*c + c_a*3*crest_quad + c_c)*G, ...
                                        (c_a*3*c + c_b)*G_quad, c_a*G_cub];
    c = c_a*c^3 + c_b*c^2 + c_c*c + c_d + (c_a*3*c + c_b)*crest_quad;
end

function [c,G,Grest,expMat,id,d] = reducePolyZono(c,G,Grest,expMat, ...
                                                                id,d,nrGen)
% reduce the number of generators of a polynomial zonotope, where we
% exploit that an interval remainder is added when reducing with Girards
% method

    n = length(c);

    if ~isempty(nrGen) && nrGen < size(G,2) + size(Grest,2)
       temp = reduce(polyZonotope(c,G,Grest,expMat,id),'girard',nrGen/n);
       c = temp.c; G = temp.G; expMat = temp.expMat; id = temp.id;
       m = max(1,size(temp.Grest,2)-n);
       d = d + rad(interval(zonotope(zeros(n,1),temp.Grest(:,m:end))));
       Grest = temp.Grest(:,1:m-1);
    end
end

function [c,G,Grest,expMat,id] = initialOrderReduction(c,G,Grest, ...
                                                       expMat,id,id_,nrGen)
% order reduction in the first layer, where we potentially redefine large
% independent generators as new dependent generators

    temp = reduce(polyZonotope(c,G,Grest,expMat,id), ...
                                                 'girard',nrGen/length(c));

    if isempty(temp.G)
        c = temp.c; G = temp.Grest; Grest = zeros(size(c));
        expMat = eye(size(G,2)); 
        id = id_ + (1:size(G,2))';
    else
        c = temp.c; G = temp.G; Grest = temp.Grest; 
        expMat = temp.expMat; id = temp.id;
    end
end

%------------- END OF CODE --------------