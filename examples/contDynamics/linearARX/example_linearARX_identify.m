function res = example_linearARX_identify
% example_linearARX_identify - example for system identification using
%   linear ARX models
%
% Syntax:
%    res = example_linearARX_identify
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean

% Authors:       Niklas Kochdumper
% Written:       30-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % trajectory data for a throttle
    y = [0 0 90 90 0 0];
    ty = [0 0.05 0.23 1 1.32 3];

    tu = [0 0.04 0.041 0.99 1 3];
    u = [0 0 0.6 0.6 0 0];

    t = 0:0.01:3;
    y = interp1(ty,y,t);
    u = interp1(tu,u,t);

    % system identification
    sys = linearARX.identify(y,t,u);

    % simulation
    simOpts.y0 = y(1:sys.n_p);
    simOpts.tFinal = t(end);
    simOpts.u = u;

    [t_,~,~,y_] = simulate(sys,simOpts);

    % visualization
    figure; hold on; box on;
    plot(t,y);
    plot(t_,y_);

    res = 1;

% ------------------------------ END OF CODE ------------------------------
