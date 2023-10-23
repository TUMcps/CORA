function example_manual_generateRandom()
% example_manual_generateRandom - example from the manual demonstrating the 
% generateRandom operation as defined in the manual
%
% Syntax:
%   example_manual_generateRandom()
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

% Authors:        Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% set seed
rng(2)

% generate random set
S = interval.generateRandom( ...
    "Dimension",2)

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:manual")
plot(S)

enlargeAxis(1.2)
title('$\mathcal{S}$','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
