function [t,x,ind,y] = simulate(obj,params,varargin)
% simulate - simulates a linear system 
%
% Syntax:
%    [t,x] = simulate(obj,params)
%    [t,x,ind] = simulate(obj,params,options)
%    [t,x,ind,y] = simulate(obj,params,options)
%
% Inputs:
%    obj - linearSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .timeStep: time step to get from initial time to final time
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
%    options - ODE45 options (optional)
%
% Outputs:
%    t - time vector
%    x - state trajectory
%    ind - returns the event which has been detected
%    y - output trajectory
%
% Example: 
%    A = [1 0; 0 2];
%    B = [1;2];
%    sys = linearSys('test',A,B);
%
%    params.x0 = [1;2];
%    params.tFinal = 2;
%
%    [t,x] = simulate(sys,params);
%
%    plot(x(:,1),x(:,2),'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       03-May-2007 
% Last update:   20-March-2008
%                08-May-2020 (MW, update interface)
%                13-February-2023 (MW, remove duplicates)                                           
% Last revision: 16-November-2021 (MW, include params.w|v)

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
isOpt = false;
if nargin >= 3 && ~isempty(varargin{1})
    options = varargin{1}; 
    isOpt = true;
end

% note: isfield-checks are required since this function can be called
% directly without any validateOptions call which would provide default
% values for the respective sets (this is only the case, if simulate is
% called via simulateRandom)

% set tStart
if ~isfield(params,'tStart')
    params.tStart = 0;
end

% shift time horizon to start at 0
tFinal = params.tFinal - params.tStart;

% set time step
if isfield(params,'timeStep')
    % time vector
    tSpan = 0:params.timeStep:tFinal;
    % append last entry if time step does not divide time horizon exactly
    if abs(tSpan(end)-tFinal) > 1e-10
        tSpan = [tSpan,tFinal];
    end
    % number of steps
    steps = length(tSpan)-1;
    timeStepGiven = true;
else
    % time vector and number of steps
    tSpan = [0,tFinal];
    steps = 1;
    timeStepGiven = false;
end

% check values of u, w, and v
[params,steps,tSpan] = aux_uwv(obj,params,steps,tSpan);
tSpan = tSpan(1:2);

% initializations
params_ = params;
t = [];
x = [];
ind = [];
y = [];
x0 = params.x0;

% computation of output set desired / possible
comp_y = nargout == 4 && ~isempty(obj.C);

% loop over all time steps
for i = 1:steps

    % input and disturbance are processed in getfcn
    params_.u = params.u(:,i);
    params_.w = params.w(:,i);

    % simulate using MATLABs ode45 function
    try
        if isOpt
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0,options);
        else
            [t_,x_,~,~,ind] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    catch
        if isOpt
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0,options);
        else
            [t_,x_] = ode45(getfcn(obj,params_),tSpan,x0);
        end
    end

    % store the results
    if i == 1
        if timeStepGiven
            % only save first and last entry
            t = [t_(1); t_(end)] + params.tStart;
            x = [x_(1,:); x_(end,:)];
        else
            t = t_ + params.tStart;
            x = x_;
        end
    else
        if timeStepGiven
            % only save last entry
            t = [t; t_(end) + t(end)];
            x = [x; x_(end,:)];
        else
            % save all simulated points
            t = [t; t_(2:end) + t(end)];
            x = [x; x_(2:end,:)];
        end
    end
    x0 = x(end,:)';

    % compute output: skip last entry since the input is not valid there,
    % instead, this will be covered by the next input
    if comp_y
        if timeStepGiven
            y_ = obj.C * x_(1,:)' + obj.D * params.u(:,i) + obj.k + params.v(:,i);
        else
            y_ = obj.C * x_(1:end-1,:)' + obj.D * params.u(:,i) + obj.k + params.v(:,i);
        end
        y = [y;y_'];
    end

    if ~isempty(ind)
        return; 
    end

end

% compute last output
if comp_y
    ylast = obj.C * x_(end,:)' + obj.D * params.u(:,end) + obj.k + params.v(:,end);
    y = [y; ylast'];
end

end


% Auxiliary functions -----------------------------------------------------

function [params,steps,tSpan] = aux_uwv(obj,params,steps,tSpan)
% set input vector u, disturbance w, and sensor noise v correctly

% set default values
if ~isfield(params,'u')
    params.u = zeros(obj.nrOfInputs,1); 
end
if ~isfield(params,'w')
    params.w = zeros(obj.dim,1);
end
if ~isfield(params,'v')
    params.v = zeros(obj.nrOfOutputs,1);
end

% check sizes
sizeU = size(params.u,2);
sizeW = size(params.w,2);
sizeV = size(params.v,2);

% if a timeStep is given, then u, w, and v have to have
%   either 1 column (fixed value) or
%   u: s (if feedthrough matrix is all-zero)
%   u: s+1 (if feedthrough matrix is non-zero)
%   w: s
%   v: s+1
% using s = number of steps

isD = any(any(obj.D));

if steps > 1 % equal to isfield(params,'timeStep')
    
    if isD
        columns_u = steps+1;
    else
        columns_u = steps;
    end

    % params.u|w|v have to match timeStep
	if sizeU ~= 1 && sizeU ~= columns_u
        throw(CORAerror('CORA:wrongFieldValue','params.u',...
            'has to be of size 1 or of the size of steps'));
    elseif sizeW ~= 1 && sizeW ~= steps
        throw(CORAerror('CORA:wrongFieldValue','params.w',...
            'has to be of size 1 or of the size of steps'));
    elseif sizeV ~= 1 && sizeV ~= steps+1
        throw(CORAerror('CORA:wrongFieldValue','params.v',...
            'has to be of size 1 or of the size of steps'));
	end
    
    % repeat fixed values
    if sizeU == 1
        params.u = repmat(params.u,1,steps+1);
    elseif sizeU == steps
        % extend by dummy value
        params.u = [params.u, zeros(size(params.u,1),1)];
    end
    if sizeW == 1
        params.w = repmat(params.w,1,steps);
    end
    if sizeV == 1
        params.v = repmat(params.v,1,steps+1);
    end


else % steps == 1, i.e., no params.timeStep given

    % check if sizes are either 1 or match one another
    
    if all([sizeU, sizeW, sizeV] == 1)
        % all sizes are 1, equals number of steps -> return
        return;
    else
        % check whether sizes fit one another
        sizesMatch = true;
        % w always one less than v
        if (sizeW ~= 1 && sizeV ~= 1) && sizeW ~= sizeV - 1
            sizesMatch = false;
        elseif isD
            if (sizeU ~= 1 && sizeV ~= 1) && sizeU ~= sizeV
                sizesMatch = false;
            elseif (sizeU ~= 1 && sizeW ~= 1) && sizeU ~= sizeW + 1
                sizesMatch = false;
            end
        else
            if (sizeU ~= 1 && sizeV ~= 1) && sizeU ~= sizeV - 1
                sizesMatch = false;
            elseif (sizeU ~= 1 && sizeW ~= 1) && sizeU ~= sizeW
                sizesMatch = false;
            end
        end
        if ~sizesMatch
            throw(CORAerror('CORA:specialError',...
                sprintf([...
                'If params.timeStep is not provided, the number of columns of the\n' ...
                '    input signal (params.u)\n' ...
                '    disturbance (params.w)\n' ...
                '    sensor noise (params.v)\n' ...
                'have to be either 1 or match one another, that is\n'...
                '    Case #1: feedthrough matrix D provided:\n'...
                '        params.u: x+1 columns\n'...
                '        params.w: x   columns\n'...
                '        params.v: x+1 columns\n'...
                '    Case #2: feedthrough matrix D not provided:\n'...
                '        params.u: x   columns\n'...
                '        params.w: x   columns\n'...
                '        params.v: x+1 columns'...
                ])));
        end
    end

    % compute steps
    if sizeV ~= 1
        steps = sizeV - 1;
    elseif sizeW ~= 1
        steps = sizeW;
    elseif sizeU ~= 1 && isD
        steps = sizeU - 1;
    elseif sizeU ~= 1 && ~isD
        steps = sizeU;
    end
    
    % lengthen vectors
    if sizeU == 1
        params.u = repmat(params.u,1,steps+1);
    elseif sizeU == steps
        % extend by dummy value
        params.u = [params.u, zeros(size(params.u,1),1)];
    end
    if sizeW == 1
        params.w = repmat(params.w,1,steps);
    end
    if sizeV == 1
        params.v = repmat(params.v,1,steps+1);
    end

    % recompute tSpan
    tSpan = linspace(0,tSpan(end),steps+1);
    tSpan = tSpan(1:2);
end

end

% ------------------------------ END OF CODE ------------------------------
