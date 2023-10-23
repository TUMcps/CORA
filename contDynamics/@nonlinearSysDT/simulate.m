function [t,x,ind,y] = simulate(obj,params)
% simulate - simulates a nonlinear discrete-time system
%
% Syntax:
%    [t,x] = simulate(obj,params)
%    [t,x,ind,y] = simulate(obj,params)
%
% Inputs:
%    obj - nonlinearSysDT object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .u: piecewise constant input signal u(t) specified as a matrix
%           for which the number of rows is identical to the number of
%           system inputs
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - returns the event which has been detected
%    y - output vector
%
% Example: 
%    f = @(x,u) [x(1) + u(1);x(2) + u(2)*cos(x(1));x(3) + u(2)*sin(x(1))];
%    dt = 0.25;
%    sys = nonlinearSysDT(f,dt);
%    
%    params.x0 = [0;0;0];
%    params.tFinal = 1;
%    params.u = [0.02 0.02 0.02 0.02; 2 3 4 5;];
%
%    [t,x] = simulate(sys,params);
%
%    plot(x(:,2),x(:,3),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT/simulate

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       22-August-2012
% Last update:   29-January-2018 (NK)
%                24-March-2020 (NK)
%                08-May-2020 (MW, update interface)
%                25-March-2021 (MA, initial state and time removed)
%                22-June-2023 (LL, keep initial state and time, add output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargout == 3
    ind = [];
end
comp_y = false;
if nargout == 4
    ind = [];
    %warning("Output trajectories not supported for class nonlinearSysDT!");
    if ~isempty(obj.out_mFile)
        comp_y = true;
    else
    	y = [];
    end
end

% parse input arguments
if ~isfield(params,'tStart')
	params.tStart = 0; 
end

t = params.tStart:obj.dt:params.tFinal;
t = t';

if t(end) ~= params.tFinal
	throw(CORAerror('CORA:specialError',...
        'Final time has to be a multiple of the sampling time!')); 
end

% consider changing inputs
change = false;

if size(params.u,2) ~= 1
    change = true;
    if (~comp_y && size(params.u,2) ~= length(t)-1) || ...
            (comp_y && size(params.u,2) ~= length(t))
        throw(CORAerror('CORA:specialError',...
            'Input signal "params.u" has the wrong dimension!'));
    end
end

% initialization
x = zeros(length(t),length(params.x0));
x(1,:) = params.x0';

if comp_y
    y = zeros(length(t),obj.nrOfOutputs);
    y(1,:) = obj.out_mFile(x(1,:)',params.u(:,1))';
end

% loop over all time steps
for i = 1:length(t)-1

    if change
        temp = obj.mFile(x(i,:)',params.u(:,i));
        if comp_y
            y(i+1,:) = obj.out_mFile(temp,params.u(:,i+1))';
        end
    else
        temp = obj.mFile(x(i,:)',params.u);
        if comp_y
            y(i+1,:) = obj.out_mFile(temp,params.u)';
        end
    end

    x(i+1,:) = temp';
end

% ------------------------------ END OF CODE ------------------------------
