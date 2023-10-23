function example_manual_specification()
% example_manual_specification - example from the manual demontrating the
% specification constructor as defined in the manual
%
% Syntax:
%   example_manual_specification()
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

% first specification
S = ellipsoid(diag([4,4]));
spec1 = specification(S,'safeSet');

% second specification
S = interval([1;1],[2.5;2.5]);
spec2 = specification(S,'unsafeSet');

% combination of both specifications
spec = add(spec1,spec2);

% plot --------------------------------------------------------------------

figure; hold on;
plot(spec)

enlargeAxis(1.2)
% title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
