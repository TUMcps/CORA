function example_manual_conHyperplane()
% example_manual_conHyperplane - example from the manual demonstrating the 
% conHyperplane constructor as defined in the manual
%
% Syntax:
%   example_manual_conHyperplane()
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

% construct constrained hyperplane
c = [1 1];
d = 1;
A = [0 1];
b = 1;

ch = conHyperplane(c,d,A,b);

% plot --------------------------------------------------------------------

figure; hold on
xlim([-1 2]); ylim([-1 2]);
plot(ch);

title('$\mathcal{CH}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
