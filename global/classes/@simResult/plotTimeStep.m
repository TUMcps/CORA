function h = plotTimeStep(obj,varargin)
% plotTimeStep - plots the time step size used in individual
%    simulations over time (all simResult simulations in one graph)
%
% Syntax:  
%    h = plotTimeStep(obj)
%    h = plotTimeStep(obj,'k',...)
%
% Inputs:
%    obj - simResult object
%    varargin - plotting preferences
%
% Outputs:
%    h - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author:       Mark Wetzlinger
% Written:      18-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% input processing
type{1} = 'b';
% If only more than one argument is passed (plotting preferences)
if nargin > 1
    type = varargin;
end

% loop over all simulations
nrSim = length(obj.t);
hold on; box on;
for i=1:nrSim
    % time axis
    cumsumtVec = [obj.t{i}(1);repelem(obj.t{i}(2:end-1),2);obj.t{i}(end)];
    % time step sizes
    tVec = repelem(diff(obj.t{i}),2);
    % plot
    plot(cumsumtVec,tVec,type{:});
end
% labels
title('Simulation: Time Step Size');
xlabel('t');
ylabel('Time Step Size');

h = get(groot,'CurrentFigure');

%------------- END OF CODE --------------