function example_manual_linComb()
% example_manual_linComb - example from the manual demonstrating the linear
%                          combination operation
%
% Syntax:  
%   example_manual_linComb()
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_manual_convHull, example_manual_enclose
%
% Author:        Niklas Kochdumper
% Written:       18-May-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % set S1 and S2
    S1 = polyZonotope([0.5;0.5],...
                      [1 1;-1 1],...
                      [],[1 2]);
    S2 = zonotope([-1.5;-1.5],...
                  [1 0;0 1]);

    % linear combination
    res = linComb(S1,S2);

    % visualization
    figure; hold on; box on;
    plot(S1,[1,2],'r');
    plot(S2,[1,2],'Color',[0 .5 0]);
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$\mathcal{S}_1 \ and \ \mathcal{S}_1$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

    figure; hold on; box on;
    plot(res,[1,2],'b','Splits',25);
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$linComb(\mathcal{S}_1,\mathcal{S}_2)$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

%------------- END OF CODE --------------