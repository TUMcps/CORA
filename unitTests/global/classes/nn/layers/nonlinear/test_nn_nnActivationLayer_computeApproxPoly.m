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
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       17-February-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

rng(1)

% test regression for all activation functions
for act = ["relu","sigmoid","tanh"]
    layer = nnActivationLayer.instantiateFromString(act);

    for poly_method =["regression", "ridgeregression"]
        for order = 1:3
            for n=1:10
                l = -3 + rand(1) * 6;
                u = l + rand(1);
        
                [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
        
                x = linspace(l, u);
                y = layer.f(x);
                y_p = polyval(coeffs, x);
    
                res = res && all(contains(interval(y_p-d, y_p+d)', y'));
            end
        end
    end
end

% test custom poly_methods
for act = ["sigmoid","tanh"]
    layer = nnActivationLayer.instantiateFromString(act);

    for poly_method =["taylor", "throw-catch"]
        for order = 1:3
            for n=1:10
                l = -3 + rand(1) * 6;
                u = l + rand(1);
        
                [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
        
                x = linspace(l, u);
                y = layer.f(x);
                y_p = polyval(coeffs, x);
    
                res = res && all(contains(interval(y_p-d, y_p+d)', y'));
            end
        end
    end
end

% test singh
for act = ["relu","sigmoid","tanh"]
    layer = nnActivationLayer.instantiateFromString(act);

    for poly_method =["singh"]
        if strcmp(act, "relu")
            orders = 1:2;
        else
            orders = 1;
        end
        for order = orders
            for n=1:2
                l = -3 + rand(1) * 6;
                u = l + rand(1);
        
                [coeffs, d] = layer.computeApproxPoly(l, u, order, poly_method);
        
                x = linspace(l, u);
                y = layer.f(x);
                y_p = polyval(coeffs, x);
    
                res = res && all(contains(interval(y_p-d, y_p+d)', y'));
            end
        end
    end
end


% ------------------------------ END OF CODE ------------------------------
