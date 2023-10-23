function example_manual_halfspace()
% example_manual_halfspace - example from the manual demonstrating the 
% halfspace constructor as defined in the manual
%
% Syntax:
%   example_manual_halfspace()
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

% construct halfspace
c = [1 1];
d = 1;

hs = halfspace(c,d);

% plot --------------------------------------------------------------------

figure; hold on
xlim([-1 2]); ylim([-1 2]);
plot(hs);

title('$\mathcal{HS}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
