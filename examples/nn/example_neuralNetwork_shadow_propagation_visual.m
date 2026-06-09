function res = example_neuralNetwork_shadow_propagation_visual()
% example_neuralNetwork_shadow_propagation_visual - visualize the 
%   propagation of a constrained zonotope through a neural networks as 
%   shadows of the same hypercube. Figure is used in [1].
%
% Syntax:
%    res = example_neuralNetwork_shadow_propagation_visual()
%
% Inputs:
%    -
%
% Outputs:
%    res - 
%
% References:
%    [1] Koller, L. "Out of the Shadows: Exploring a Latent Space for 
%       Neural Network Verification". (2025) arXiv
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       21-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

rng('default')

% Specify the input set.
c = [0; 0];
G = 1/sqrt(2)*[1 0 0; 0 -1 0];
X = zonotope(c,G);

% Specify the offset for zonotope.
offset = -3;

% Specify affine map (90° rotation).
W = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
b = [1; 0];

% Specify image enclosure.
m = [1; 1/2];
dl = [0; 1/2];
du = [0; 0];

% Propagate sets.
H0 = X; % input set
H1 = W*H0 + b; % affine map
H2 = diag(m)*H1 + 1/2*(du + dl);
H2.G(2,3) = 1/2*(du(2) - dl(2)); % image enclosure

% Compute the true output set.
Y = H1 & polytope(-eye(2),zeros(2,1));

% Specify output constraints.
A = [-1 0];
b = -1.5;
apprErr = sum(abs(A*H2.G(:,3:end)));

% Falsify: Sample a potential counter example.
[~,q] = size(H2.G);
factors = sign(-A*H2.G)';

% Refine: Convert constraints to constraint on the zonotope.
C = -A*H2.G;
d = -(b - A*H2.c + apprErr);
% Add ReLU split constraint.
% C = [C; -H1.G(2,:)];
% d = [d; 0];
% Add ReLU tightening constraint.
% C = [C; -H2.G(2,:); H1.G(2,:) - H2.G(2,:)];
% d = [d; H2.c(2,:); H2.c(2,:) - H1.c(2,:)];

% Plot everything in 3D.

figure; 
subplot(1,3,1); 
hold on; box on; grid on;
title('Input Space')
xlim([-1 1])
ylim([-1 1])
zlim([offset 1])
Ti = aux_plot3DZonotope(H0,C,d,[],[],eye(3),apprErr,factors, ...
    offset,CORAcolor('CORA:reachSet'),'Zonotope');
% legend

subplot(1,3,2); 
hold on; box on; grid on;
title('Hidden Space')
xlim([-0.25 2.25])
ylim([-1.25 1.25])
zlim([offset 1])
Ti = aux_plot3DZonotope(H1,C,d,[],[],Ti,apprErr,factors, ...
    offset,CORAcolor('CORA:reachSet'),'Zonotope');
% legend

subplot(1,3,3); 
hold on; box on; grid on;
title('Output Space')
xlim([-0.25 2.25])
ylim([-1 1.5])
zlim([offset 1])
Ti = aux_plot3DZonotope(H2,C,d,A,b,Ti,apprErr,factors, ...
    offset,CORAcolor('CORA:reachSet'),'Zonotope');
% legend

% Set the output.
res = true;

end


% Auxiliary functions -----------------------------------------------------

function Ti = aux_plot3DZonotope(Z,C,d,A,b,T0,apprErr, ...
    factors,offset,color,zname)
    % Obtain number of generators.
    [~,q] = size(Z.G);

    % Extract the center.
    c = Z.c;
    
    % Construct generator space.
    B = zonotope(zeros(3,1),blkdiag(eye(q),zeros(3-q)));
    
    % Compute the transformation matrix for the hypercube.
    options = optimoptions('fsolve', ...
        Algorithm='levenberg-marquardt',Display='off');
    Ti = fsolve(@(T) reshape(eye(2,3)*T - Z.G,[],1),T0,options);
    % Scale last dimension.
    % Ti(3,:) = Ti(3,:)*1/2*sqrt(2);

    % Transform the unit hypercube s.t. its shadow is the given zonotope.
    B = [c; 0] + Ti*B;
    % Convert constraint to equality constraints.
    [G,C,d] = aux_inequToEquConstr(B.G,C,d);
    % Convert constraints to equality constraints by adding slack variables.
    cB = conZonotope(B.c,G,C,d);

    % Compute the coordinates of the sample.
    x = B.c + B.G*factors;

    % Compute all vertices of the generator hypercube.
    vsZ = aux_computeProjectedVerticsOfHypercube(B);

    % Plot dashed lines between the vertices of the hypercube and the 
    % vertices of the zonotope.
    for i=1:size(vsZ,2)
        plot3([vsZ(1,i) vsZ(1,i)],[vsZ(2,i) vsZ(2,i)],[offset vsZ(3,i)], ...
            ':k','HandleVisibility','off');
    end

    % Plot the hypercube.
    plot(B,1:3,'--', ...
        'FaceColor',CORAcolor('CORA:simulations'),'FaceAlpha',0.0, ...
        'EdgeColor',CORAcolor('CORA:simulations'), ...
        'DisplayName','Hypercube' ...
    );
    % Plot the constraint hypercube.
    plot(cB,1:3, ...
        'FaceColor',CORAcolor('CORA:simulations'),'FaceAlpha',0.2, ...
        'EdgeColor',CORAcolor('CORA:simulations'),'LineWidth',2, ...
        'DisplayName','constraint Hypercube' ...
    );

    % Plot the zonotope.
    plot(B,1:2,'--', ...
        'FaceColor',color,'FaceAlpha',0.0,...
        'EdgeColor',color,'LineWidth',1, ...
        'ZPos',offset,'DisplayName',zname ...
    );
    % Plot the constraint zonotope.
    plot(cB,1:2, ...
        'FaceColor',color,'FaceAlpha',0.2,...
        'EdgeColor',color,'LineWidth',2, ...
        'ZPos',offset,'DisplayName',['constraint ' zname] ...
    );

    % Plot the unsafe set specification.
    if ~isempty(A)
        % Construct the halfspace specification.
        spec = polytope(A,b);
        % Add approximation error.
        specApprErr = polytope(A,b + apprErr);
    
        % Plot the specification.
        plot(and(spec,interval([-10;-10],[10; 10])),1:2, ...
            'FaceColor',CORAcolor('CORA:unsafe'),'FaceAlpha',0.2, ...
            'EdgeColor',CORAcolor('CORA:unsafe'),'LineWidth',2, ...
            'ZPos',offset,'DisplayName','Specification');
        % Plot the specification with approximation error.
        plot(and(specApprErr,interval([-10;-10],[10; 10])),1:2,'--', ...
            'FaceColor',CORAcolor('CORA:unsafe'),'FaceAlpha',0.0, ...
            'EdgeColor',CORAcolor('CORA:unsafe'),'LineWidth',1, ...
            'ZPos',offset,'DisplayName','Approximation error');
    end
    
    % Plot the sample in the hypercube.
    scatter3([x(1) x(1)],[x(2) x(2)],[x(3) offset],100,'.k', ...
        'DisplayName','Sample');
    % Plot the samples in the zonotope.
    plot3([x(1) x(1)],[x(2) x(2)],[offset x(3)],'--k', ...
        'HandleVisibility','off');

    % Plot chain of 3D generators.
    x_ = B.c;
    for i=1:q
        % Extract and scale the generator.
        gi = factors(i)*B.G(:,i);
        % Plot the generator.
        quiver3(x_(1),x_(2),x_(3),gi(1),gi(2),gi(3), ...
            'Color',CORAcolor('CORA:simulations'),'linewidth',1.5, ...
            'AutoScale','off','HandleVisibility','off');
        quiver3(x_(1),x_(2),offset,gi(1),gi(2),0, ...
            'Color',CORAcolor('CORA:simulations'),'linewidth',1.5, ...
            'AutoScale','off','HandleVisibility','off');
        % Update current point.
        x_ = x_ + gi;
    end
end

function vsZ = aux_computeProjectedVerticsOfHypercube(B)
    % Obtain number of dimensions.
    n = length(B.c);
    % Compute all vertices of the generator hypercube.
    vsB = vertices(B);
    % Compute the vertices of the zonotope by project the vertices of 
    % the hypercube onto the x-y-plane and only keep vertices from the
    % x-y-plane.
    idx = convhull(vsB(1:2,:)');
    vsZ = vsB(:,idx(1:(end-1)));
    % Sort the vertices for comparison.
    [~,idx] = sort(10.^(0:(n-1))*vsZ);
    vsZ = vsZ(:,idx);
end

function [G,C,d] = aux_inequToEquConstr(G,C,d)
    % Obtain number of dimensions.
    [n,~] = size(G);
    % Obtain number of constrains.
    [p,~] = size(C);
    % Convert constraints to equality constraints by adding slack 
    % variables.
    G = [G zeros([n p])];
    % Compute scale for the slack variable.
    s = 1/2*(sum(abs(C),2) + d);
    C = [C eye(p).*s];
    % Compensate for the slack variable.
    d = d - s;
end

% ------------------------------ END OF CODE ------------------------------
