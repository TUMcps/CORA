function completed = example_polyZonotope_relaxExponents()
% example_polyZonotope_relaxExponents - shows the paper example
%
% Syntax:
%    completed = example_polyZonotope_relaxExponents()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       25-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init polyZonotope
pZ = polyZonotope([2;2],[1 -1; 2 3],[],[1 1; 4 2]);

% convert to zonotope
Z = zonotope(pZ);

% relax eponents
pZ_relax = pZ.relaxExponents(1);

% convert relaxed to zonotope
Z_relax = zonotope(pZ.relaxExponents);

% plot
figure; hold on;
plot(pZ,1:2,'DisplayName','$\mathcal{PZ}$')
plot(Z,1:2,'DisplayName','$\texttt{zonotope}(\mathcal{PZ})$')
plot(pZ_relax,1:2,'DisplayName','$\texttt{relax}(\mathcal{PZ})$')
plot(Z_relax,1:2,'--','DisplayName','$\texttt{zonotope}(\texttt{relax}(\mathcal{PZ}))$')
legend('Interpreter','latex')

completed = true;

end

% ------------------------------ END OF CODE ------------------------------
