function example_manual_cartProd()
% example_manual_mtimes - example from the manual demonstrating the Cartesian Product
% sum of two sets.
%
% Syntax:  
%   example_manual_cartProd()
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
    S1 = interval(-2,1);
    S2 = interval(-1,2);
    
    % Cartesian product
    res = cartProd(S1,S2);
    
    figure; hold on; box on;
    plot(res,[1,2],'b');
    axis equal;
    xlim([-3,3]);ylim([-3,3]);
    xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
    title('$\mathcal{S}_1 \times \mathcal{S}_2$','Interpreter','latex');
    xticks(-3:3); yticks(-3:3);

%------------- END OF CODE --------------