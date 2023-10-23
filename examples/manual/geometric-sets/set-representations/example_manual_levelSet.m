function example_manual_levelSet()
% example_manual_levelSet - example from the manual demonstrating the 
% levelSet constructor as defined in the manual
%
% Syntax:
%   example_manual_levelSet()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct level set
vars = sym('x',[2,1]);
f = 1/vars(1)^2 - vars(2);
op = '==';

ls = levelSet(f,vars,op);

% plot --------------------------------------------------------------------

figure; hold on
xlim([0.5 3]); ylim([0 3]);
plot(ls);
xlim([0 3]);

title('$\mathcal{LS}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
