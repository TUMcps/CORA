function example_manual_mtimes()
% example_manual_mtimes - example from the manual demonstrating the linear
% map of a set as defined in the manual.
%
% Syntax:  
%   example_manual_mtimes()
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

    % set S
    S = zonotope([0 1 1 0;...
        0 1 0 1]);
    M = [1 0; -1 0.5];
    
    % mtimes
    res = M * S;
    
    figure; hold on; box on;
    plot(S,[1,2],'r');
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$\mathcal{S}$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);
    
    figure; hold on; box on;
    plot(res,[1,2],'b');
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$M \otimes \mathcal{S}$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

%------------- END OF CODE --------------