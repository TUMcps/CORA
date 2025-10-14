function [t,x,ind,y] = simulate(linsysDT,params,varargin)
% simulate - simulates a linear discrete-time system
%
% Syntax:
%    [t,x] = simulate(linsysDT,params)
%    [t,x,ind,y] = simulate(linsysDT,params)
%
% Inputs:
%    linsysDT - linearSysDT object
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
%       .w: disturbance of dimension r x s (r = disturbance dimension)
%       .v: sensor noise of dimension p x s+1 (p = noise dimension)
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
%    linsysDT = linearSysDT(A,B,dt);
%
%    params.x0 = [0;0];
%    params.tFinal = 4;
%    params.u = [0.1 0.05 0.05 -0.1];
%
%    [t,x] = simulate(linsysDT,params);
%
%    plot(x(1,:),x(2,:),'.k','MarkerSize',20);
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
%                28-August-2025 (LL, transpose t, x, and y)
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
t = params.tStart:linsysDT.dt:params.tFinal;
steps = length(t)-1;

% check end of time vector
if t(end) ~= params.tFinal
    throw(CORAerror('CORA:specialError',...
        'Final time has to be a multiple of the sampling time!'));
end

% check u, w, and v; extend fixed vector to matrices for loop below
params = aux_uwv(linsysDT,params,steps);

% initial state
x = zeros(linsysDT.nrOfDims,steps+1);
x(:, 1) = params.x0;

% computation of output set desired / possible
comp_y = nargout == 4 && ~isempty(linsysDT.C);

if ~iscell(linsysDT.A)
    % initial output
    if comp_y
    	y = zeros(linsysDT.nrOfOutputs,steps+1);
        y(:,1) = linsysDT.C * x(:,1) + linsysDT.D * params.u(:,1) + linsysDT.k + linsysDT.F * params.v(:,1);
    end

    % loop over all time steps
    for i = 1:steps

        % compute successor state
        temp = linsysDT.A * x(:,i) + linsysDT.B * params.u(:,i) + linsysDT.c + linsysDT.E * params.w(:,i);
        x(:,i+1) = temp';

        % compute output
        if comp_y
            y(:,i+1) = linsysDT.C * x(:,i+1) + linsysDT.D * params.u(:,i+1) + linsysDT.k + linsysDT.F * params.v(:,i+1);
        end
    end
else
    if comp_y
    	y = zeros(linsysDT.nrOfOutputs,steps+1);
        y(:,1) = linsysDT.C{1} * x(:,1) + linsysDT.D{1} * params.u(:,1) + linsysDT.k + linsysDT.F * params.v(:,1);
    end

    % loop over all time steps
    for i = 1:steps

        % compute successor state
        temp = linsysDT.A{i} * x(:,i) + linsysDT.B{i} * params.u(:,i) + linsysDT.c + linsysDT.E * params.w(:,i);
        x(:,i+1) = temp';

        % compute output
        if comp_y
            y(:,i+1) = linsysDT.C{i+1} * x(:,i+1) + linsysDT.D{i+1} * params.u(:,i+1) + linsysDT.k + linsysDT.F * params.v(:,i+1);
        end
    end
end
end


% Auxiliary functions -----------------------------------------------------

function params = aux_uwv(sys,params,steps)

% check input, set default
if ~isfield(params,'u')
    % should be of length steps if no output equation is given and steps+1
    % if an output equation is given, but we keep it simple
    params.u = zeros(sys.nrOfInputs,steps+1);
else
    % length depends on whether feedthrough matrix is all-zero or not
    if (~iscell(sys.D) && any(any(sys.D))) || (iscell(sys.D))
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
    elseif (~iscell(sys.D) && ~any(any(sys.D)))
        % extend by dummy value; note that this value should not have any
        % effects on the result as D is all-zero
        params.u = [params.u, zeros(size(params.u,1),1)];
    end
end

% check disturbance, set default
if ~isfield(params,'w')
    params.w = zeros(sys.nrOfDisturbances,steps);
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
    params.v = zeros(sys.nrOfNoises,steps+1);
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
