function res = enlargeAxis(factor)
% enlargeAxis - enlarges the axis of the current figure
%
% Syntax:
%    enlargeAxis()
%    enlargeAxis(factor)
%
% Inputs:
%    factor - enlargement factor
%
% Outputs:
%    res
%

% Authors:       Tobias Ladner
% Written:       09-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 1
    factor = 1.2;
end

% get axis limits
xLim = xlim();
yLim = ylim();

% enlarge viewbox
I = interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]);
I = enlarge(I,factor);

% set axis
xlim([I.inf(1) I.sup(1)])
ylim([I.inf(2) I.sup(2)])


end

% ------------------------------ END OF CODE ------------------------------
