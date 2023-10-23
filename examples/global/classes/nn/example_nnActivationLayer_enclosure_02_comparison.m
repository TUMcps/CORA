function res = example_nnActivationLayer_enclosure_02_comparison()
% example_nnActivationLayer_enclosure_02_comparison - example for 
%    neural network enclosure of an activation layer using different
%    polynomial methods
%
% Syntax:
%    res = example_nnActivationLayer_enclosure_02_comparison()
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
methods = ["regression", "singh"];

l_max = l-1;
u_max = u+1;

figure;
sgtitle("Polynomial Methods Comparison")

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
    useCORAcolors("CORA:default")

    x = linspace(l, u);

    for j=1:length(methods)
        poly_method = methods(j);
        % compute enclosure
        [coeffs, d] = layer.computeApproxPoly(l, u, 1, poly_method);
        y_p = polyval(coeffs, x);

        % plot incl. bounds
        color = CORAcolor('CORA:next');
        plot(x, y_p, "Color", color, 'DisplayName', poly_method);
        plot(x, y_p-d, '--', "Color", color, 'HandleVisibility', 'off');
        plot(x, y_p+d, '--', "Color", color, 'HandleVisibility', 'off');
        updateColorIndex(j)
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
