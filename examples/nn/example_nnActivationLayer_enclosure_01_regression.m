function res = example_nnActivationLayer_enclosure_01_regression()
% example_nnActivationLayer_enclosure_01_regression - example for 
%    neural network enclosure of an activation layer using regression
%
% Syntax:
%    res = example_nnActivationLayer_enclosure_01_regression()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
l = -1;
u = 2;
max_order = 5;

l_max = l-1;
u_max = u+1;
cs = hsv(max_order);


figure;
sgtitle("Enclosure using Regression")

acts = ["Relu", "Sigmoid", "Tanh"];
for i=1:length(acts)
    act = acts(i);
    layer = nnActivationLayer.instantiateFromString(act);

    subplot(1, length(acts), i); hold on;
    title(act);

    % plot f(x)
    x_max = linspace(l_max, u_max);
    y_max = layer.f(x_max);
    plot(x_max, y_max, 'k', 'DisplayName', "f(x)")

    x = linspace(l, u);

    for order=1:max_order
        % compute enclosure
        [coeffs, d] = layer.computeApproxPoly(l, u, order, 'regression');
        y_p = polyval(coeffs, x);

        % plot incl. bounds
        plot(x, y_p, "Color", cs(order, :), 'DisplayName', sprintf("Order: %d", order));
        plot(x, y_p-d, '--', "Color", cs(order, :), 'HandleVisibility', 'off');
        plot(x, y_p+d, '--', "Color", cs(order, :), 'HandleVisibility', 'off');
    end

    % plot [l, u] limits
    gray = ones(1, 3) * 0.7;
    plot([l,l],ylim, '--', "Color", gray, 'HandleVisibility', 'off')
    plot([u,u],ylim, '--', "Color", gray, 'HandleVisibility', 'off')

    % plot legends
    legend('Location', 'best')
end


res = true;

end

% ------------------------------ END OF CODE ------------------------------
