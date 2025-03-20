function res = test_nn_nnActivationLayer_computeApproxPoly()
% test_nn_nnActivationLayer_computeApproxPoly - tests the polynomial
%    approximation of the activation function
%
% Syntax:
%    res = test_nn_nnActivationLayer_computeApproxPoly()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       17-February-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test regression for all activation functions
for act = ["relu","sigmoid","tanh"]
    layer = nnActivationLayer.instantiateFromString(act);

    for poly_method =["regression", "ridgeregression"]
        for order = 1:3
            for n=1:10
                % init bounds
                l = -3 + rand(1) * 6;
                u = l + rand(1);
        
                % compute coefficients and error
                [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
        
                % evaluate polynomial
                x = linspace(l, u);
                y = layer.f(x);
                y_p = polyval(coeffs, x);
    
                % assert containment
                assertLoop(all(contains(interval(y_p-d, y_p+d)', y')),'',[],act,poly_method,order,n);
            end
        end
    end
end

% test custom poly_methods
for act = ["sigmoid", "tanh"]
    layer = nnActivationLayer.instantiateFromString(act);

    for poly_method =["taylor", "throw-catch"]
        for order = 1:3
            for n=1:10
                % init bounds
                l = -3 + rand(1) * 6;
                u = l + rand(1);
        
                % compute coefficients
                [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
        
                % evaluate polynomial
                x = linspace(l, u);
                y = layer.f(x);
                y_p = polyval(coeffs, x);
    
                % assert containment
                assertLoop(all(contains(interval(y_p-d, y_p+d)', y')),'',[],act,poly_method,order,n);
            end
        end
    end
end

% test singh
for act = ["relu","sigmoid","tanh"]
    layer = nnActivationLayer.instantiateFromString(act);
    poly_method = "singh";

    if strcmp(act, "relu")
        orders = 1:2;
    else
        orders = 1;
    end
    for order = orders
        for n=1:2
            % init bounds
            l = -3 + rand(1) * 6;
            u = l + rand(1);
    
            % compute coefficients
            [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
    
            % evaluate polynomial
            x = linspace(l, u);
            y = layer.f(x);
            y_p = polyval(coeffs, x);

            % assert containment
            assertLoop(all(contains(interval(y_p-d, y_p+d)', y')),'',[],act,poly_method,order,n);
        end
    end
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
