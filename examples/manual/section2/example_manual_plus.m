function example_manual_plus()
% example_manual_mtimes - example from the manual demonstrating the Minkowski sum of two sets.
%
% Syntax:  
%   example_manual_plus()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:
%
% Author:        Niklas Kochdumper, Philipp Gassert
% Written:       04-June-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % set S1 and S2
    S1 = zonotope([0 0.5 1;...
        0 1   0]);
    S2 = zonotope([0 1 0;...
        0 0 1]);

    % Minkowski sum
    res = S1 + S2;

    figure; hold on; box on;
    plot(S1,[1,2],'r');
    plot(S2,[1,2],'Color', [0 0.5 0]);
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$\mathcal{S}_1 \ and \ \mathcal{S}_1$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

    figure; hold on; box on;
    plot(res,[1,2],'b');
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$\mathcal{S}_1 \oplus \mathcal{S}_1$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

%------------- END OF CODE --------------