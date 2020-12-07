function [t,x] = simulate(obj,params)
% simulate - simulates a linear discrete-time system
%
% Syntax:  
%    [t,x] = simulate(obj,params)
%
% Inputs:
%    obj - nonlinearSysDT object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .u: piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system input
%
% Outputs:
%    t - time vector
%    x - state vector
%
% Example:
%    A = [1 2; -3 1];
%    B = [2;1];
%    dt = 1;
%    sys = linearSysDT(A,B,dt);
%
%    params.x0 = [0;0];
%    params.tFinal = 4;
%    params.u = [0.1 0.05 0.05 -0.1];
%
%    [t,x] = simulate(sys,params);
%
%    plot(x(:,1),x(:,2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT/simulate

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      20-March-2020
% Last update:  24-March-2020 (NK)
%               08-May-2020 (MW, update interface)
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if ~isfield(params,'tStart')
       params.tStart = 0; 
    end

    t = params.tStart:obj.dt:params.tFinal;
    t = t';

    if t(end) ~= params.tFinal
       error('Final time has to be a multiple of the sampling time!'); 
    end

    % consider changing inputs
    change = 0;

    if size(params.u,2) ~= 1
        change = 1;
        if size(params.u,2) ~= length(t)-1
            error('Input signal "params.u" has the wrong dimension!'); 
        end
    end

    % initialization
    x = zeros(length(t),length(params.x0));
    x(1,:) = params.x0';

    % get constant offset
    if isempty(obj.c)
        c = zeros(size(obj.A,1),1); 
    else
        c = obj.c;
    end

    % loop over all time steps
    for i = 1:length(t)-1

        if change
            temp = obj.A*x(i,:)' + obj.B * params.u(:,i) + c;
        else
            temp = obj.A*x(i,:)' + obj.B * params.u + c;
        end

        x(i+1,:) = temp';
    end
    
%------------- END OF CODE --------------