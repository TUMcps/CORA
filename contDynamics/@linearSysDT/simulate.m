function [t,x,ind,y] = simulate(obj,params)
% simulate - simulates a linear discrete-time system
%
% Syntax:
%    [t,x] = simulate(obj,params)
%    [t,x,ind,y] = simulate(obj,params)
%
% Inputs:
%    obj - linearSysDT object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .u: piecewise-constant input signal u(t) of dimension
%               m x z, where
%               m = number of system inputs, and
%               z = 1 (fixed input) or
%                   s (if no feedthrough matrix given)
%                   s+1 (if feedthrough matrix given)
%                   where s is the number of steps
%       .w: disturbance of dimension n x s (n = system dimension)
%       .v: sensor noise of dimension n_out x s+1 (n_out = output dimension)
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - [] (argument exists only for syntax consistency)
%    y - output vector
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

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       20-March-2020
% Last update:   24-March-2020 (NK)
%                08-May-2020 (MW, update interface)
%                12-January-2021 (MA, disturbance input added)
% Last revision: 16-November-2021 (MW, simplify loop, rewrite handling of params.u|w|v, integrate output)

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if ~isfield(params,'tStart')
	params.tStart = 0; 
end

% set default output arguments
ind = [];
y = [];

% compute time vector and number of steps
t = (params.tStart:obj.dt:params.tFinal)';
steps = length(t)-1;

% check end of time vector
if t(end) ~= params.tFinal
    throw(CORAerror('CORA:specialError',...
        'Final time has to be a multiple of the sampling time!'));
end

% check u, w, and v; extend fixed vector to matrices for loop below
params = aux_uwv(obj,params,steps);

% initial state
x = zeros(steps+1,obj.dim);
x(1,:) = params.x0';

% computation of output set desired / possible
comp_y = nargout == 4 && ~isempty(obj.C);

% initial output
if comp_y
	y = zeros(steps+1,obj.nrOfOutputs);
    y(1,:) = obj.C * x(1,:)' + obj.D * params.u(:,1) + obj.k + params.v(:,1);
end


% loop over all time steps
for i = 1:steps
    
    % compute successor state
    temp = obj.A * x(i,:)' + obj.B * params.u(:,i) + obj.c + params.w(:,i);
    x(i+1,:) = temp';

    % compute output
    if comp_y
        y(i+1,:) = obj.C * x(i+1,:)' + obj.D * params.u(:,i+1) + obj.k + params.v(:,i+1);
    end
    
end

end


% Auxiliary functions -----------------------------------------------------

function params = aux_uwv(obj,params,steps)

% check input, set default
if ~isfield(params,'u')
    % should be of length steps if no output equation is given and steps+1
    % if an output equation is given, but we keep it simple
    params.u = zeros(obj.nrOfInputs,steps+1);
else
    % length depends on whether feedthrough matrix is all-zero or not
    if any(any(obj.D))
        columns_u = steps+1;
    else
        columns_u = steps;
    end
	if size(params.u,2) == 1
        % length always steps+1
        params.u = repmat(params.u,1,steps+1);
    elseif size(params.u,2) ~= columns_u && size(params.u,2) ~= steps+1
        throw(CORAerror('CORA:wrongFieldValue','params.u',...
            'has to be of size(2) of steps+1'));
    elseif ~any(any(obj.D))
        % extend by dummy value; note that this value should not have any
        % effects on the result as D is all-zero
        params.u = [params.u, zeros(size(params.u,1),1)];
    end
end

% check disturbance, set default
if ~isfield(params,'w')
    params.w = zeros(obj.dim,steps);
else
    if size(params.w,2) == 1
        params.w = repmat(params.w,1,steps);
    elseif size(params.w,2) ~= steps
        throw(CORAerror('CORA:wrongFieldValue','params.w',...
            'has to be of size(2) of steps+1'));
    end
end

% check sensor noise, set default
if ~isfield(params,'v')
    params.v = zeros(obj.nrOfOutputs,steps+1);
else
    if size(params.v,2) == 1
        params.v = repmat(params.v,1,steps+1);
    elseif size(params.v,2) ~= steps+1
        throw(CORAerror('CORA:wrongFieldValue','params.v',...
            'has to be of size(2) of steps'));
    end
end

end
    
% ------------------------------ END OF CODE ------------------------------
