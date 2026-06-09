function res = example_conZonotope_shadow_visualization()
% example_conZonotope_shadow_visualization - visualize a constrained
%   zonotope as the shadow of a constrained hypercube. Frigure is used in
%   [1].
%
% Syntax:
%    res = example_conZonotope_shadow_visualization()
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

% Create 3D hypercube.
Bc = zeros(3,1);
BG = 1/3*eye(3);
% Obtain number of generators.
q = size(BG,2);
% Rotate hypercube.
rx = pi/4;
ry = pi/4;
rz = pi/4;
Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
Ry = [cos(ry) 0 -sin(ry); 0 1 0; sin(ry) 0 cos(ry)];
Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
BG = Rx*Ry*Rz*BG;
% % Construct hypercube zonotope.
B = zonotope(Bc,BG);
% Construct constraints.
A = [BG(2,:) 1];
b = -0.1 + 1;

% Specify affine map.
M = diag([1; 3/4; 1])*Rx*Ry*Rz;

% Sample a point in the hypercube.
factors = 2*rand([q 1]) - 1;

% Offset for 2D zonotope.
offset = -1;

figure; 
subplot(1,3,1); hold on; box on; grid on;
title('Zonotope')
xlim([-1 1])
ylim([-1 1])
zlim([offset 1])
% Plot everything.
aux_plot2DZonotope(B,A,b,factors, ...
    offset,CORAcolor('CORA:reachSet'),'Zonotope');
legend

subplot(1,3,2); hold on; box on; grid on;
title('Zonotope (Projection)')
xlim([-1 1])
ylim([-1 1])
zlim([offset 1])
% Plot everything.
aux_plot3DZonotope(B,A,b, ...
    factors,offset,CORAcolor('CORA:reachSet'),'Zonotope');
legend

subplot(1,3,3); hold on; box on; grid on;
title('transformed Zonotope')
xlim([-1 1])
ylim([-1 1])
zlim([offset 1])
% Plot the transformed zonotope.
aux_plot3DZonotope(M*B,A,b,factors, ...
    offset,CORAcolor('CORA:highlight1'),'trans. Zonotope');
legend

% Set the output.
res = true;

end


% Auxiliary functions -----------------------------------------------------

function aux_plot2DZonotope(B,A,b,factors,offset,color,zname)
    % Compute the sampled point.
    x3 = B.c + B.G*factors;

    % Construct constraint zonotope.
    cB = conZonotope(B.c,[B.G zeros(3,1)],A,b);

    % Obtain number of generators.
    q = size(B.G,2);
    
    % Plot the 2D zonotope.
    plot(B,1:2,'--', ...
        'FaceColor',color,'FaceAlpha',0.0,...
        'EdgeColor',color,'LineWidth',1, ...
        'DisplayName',zname ...
    )
    % Plot the 2D constraint zonotope.
    plot(cB,1:2, ...
        'FaceColor',color,'FaceAlpha',0.2,...
        'EdgeColor',color,'LineWidth',2, ...
        'DisplayName',['constraint ' zname] ...
    )
    
    % Plot the sample.
    scatter3([x3(1) x3(1)],[x3(2) x3(2)],[x3(3) offset],100,'.k', ...
        'DisplayName','Sample')
    plot3([x3(1) x3(1)],[x3(2) x3(2)],[offset x3(3)],'--k', ...
        'HandleVisibility','off')
    % Plot chain of 3D generators.
    x3_ = B.c;
    for i=1:q
        % Extract and scale the generator.
        gi = factors(i)*B.G(:,i);
        % Plot the generator.
        quiver3(x3_(1),x3_(2),offset,gi(1),gi(2),0, ...
            'Color',CORAcolor("CORA:simulations"),'linewidth',1.5, ...
            'AutoScale','off','HandleVisibility','off');
        % Update current point.
        x3_ = x3_ + gi;
    end
end

function aux_plot3DZonotope(B,A,b,factors,offset,color,zname)
    % Compute the sampled point.
    x3 = B.c + B.G*factors;

    % Construct constraint zonotope.
    cB = conZonotope(B.c,[B.G zeros(3,1)],A,b);

    % Obtain number of generators.
    q = size(B.G,2);
    
    % Compute all vertices.
    vs = vertices(B);
    idx = convhull(vs(1:2,:)');
    vs = vs(:,idx);

    % Plot dashed lines between the vertices.
    for i=1:size(vs,2)
        plot3([vs(1,i) vs(1,i)],[vs(2,i) vs(2,i)],[offset vs(3,i)],':k', ...
        'HandleVisibility','off')
    end
    % Plot the hypercube.
    plot(B,1:3,'--', ...
        'FaceColor',CORAcolor('CORA:simulations'),'FaceAlpha',0.0, ...
        'EdgeColor',CORAcolor('CORA:simulations'), ...
        'DisplayName','Hypercube' ...
    )
    plot(cB,1:3, ...
        'FaceColor',CORAcolor('CORA:simulations'),'FaceAlpha',0.2, ...
        'EdgeColor',CORAcolor('CORA:simulations'),'LineWidth',2, ...
        'DisplayName','constraint Hypercube' ...
    )

    % Plot the 2D zonotope.
    plot(B,1:2,'--', ...
        'FaceColor',color,'FaceAlpha',0.0,...
        'EdgeColor',color,'LineWidth',1, ...
        'ZPos',offset,'DisplayName',zname ...
    )
    % Plot the 2D constraint zonotope.
    plot(cB,1:2, ...
        'FaceColor',color,'FaceAlpha',0.2,...
        'EdgeColor',color,'LineWidth',2, ...
        'ZPos',offset,'DisplayName',['constraint ' zname] ...
    )
    
    % Plot the sample.
    scatter3([x3(1) x3(1)],[x3(2) x3(2)],[x3(3) offset],100,'.k', ...
        'DisplayName','Sample')
    plot3([x3(1) x3(1)],[x3(2) x3(2)],[offset x3(3)],'--k', ...
        'HandleVisibility','off')
    % Plot chain of 3D generators.
    x3_ = B.c;
    for i=1:q
        % Extract and scale the generator.
        gi = factors(i)*B.G(:,i);
        % Plot the generator.
        quiver3(x3_(1),x3_(2),x3_(3),gi(1),gi(2),gi(3), ...
            'Color',CORAcolor("CORA:simulations"),'linewidth',1.5, ...
            'AutoScale','off','HandleVisibility','off');
        quiver3(x3_(1),x3_(2),offset,gi(1),gi(2),0, ...
            'Color',CORAcolor("CORA:simulations"),'linewidth',1.5, ...
            'AutoScale','off','HandleVisibility','off');
        % Update current point.
        x3_ = x3_ + gi;
    end
end

% ------------------------------ END OF CODE ------------------------------
