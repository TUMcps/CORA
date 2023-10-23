function example_manual_interval()
% example_manual_interval - example from the manual demonstrating the 
% interval constructor as defined in the manual
%
% Syntax:
%   example_manual_interval()
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

% construct interval
lb = [-2; -1];
ub = [4; 3];

int = interval(lb,ub);

% plot --------------------------------------------------------------------

figure; hold on
plot(int);

enlargeAxis(1.2)
title('$\mathcal{I}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
