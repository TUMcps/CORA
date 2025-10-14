function p = estimateParameter(sys,varargin)
% estimateParameter - Determine the values for the parameters that best fit
%   the given trajectory data
%
% Syntax:
%    p = estimateParameter(sys,traj)
%    p = estimateParameter(sys,traj,options)
%    p = estimateParameter(sys,x,t)
%    p = estimateParameter(sys,x,t,options)
%    p = estimateParameter(sys,x,t,u)
%    p = estimateParameter(sys,x,t,u,options)
%
% Inputs:
%    sys - nonlinParamSys object
%    traj - object of class "trajectory" storing the trajectory data
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,n]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [N,1]
%    u - cell-array storing the inputs for the simulated trajectories
%        as vectors of dimension [N-1,m]
%    options - algorithm options for parameter estimation
%
%       -.alg:  algorithm that is used 
%               ('singleStep' (default) or 'multiStep')
%
% Outputs:
%    p - estimated value for the system parameters
%
% Example: 
%    l = 2.8;
%    f = @(x,u,p) [x(1)*cos(x(2)); ...
%                  x(1)/p(1) * tan(u(1))];
%    sys = nonlinParamSys(f);
%
%    simOpts.x0 = [1; 0];
%    simOpts.tFinal = 1;
%    simOpts.u = [0.3];
%    simOpts.p = l;
%    [t,x] = simulate(sys,simOpts);
%    u = repmat(simOpts.u,[1,length(t)]);
%
%    p = estimateParameter(sys,x,t,u)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       12-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse the input arguments and bring trajectory data to suitable form
    [traj,options] = checkDataIdentification(varargin{:});

    if size(traj(1).x,1) ~= sys.nrOfDims
        throw(CORAerror('CORA:specialError', ...
                  'Number of states for system and data does not match!'));
    else
        m = 1;
        if isfield(traj(1),'u')
            m = size(traj(1).u,1);
        end
        if m ~= sys.nrOfInputs
            throw(CORAerror('CORA:specialError', ...
                  'Number of inputs for system and data does not match!'));
        end
    end

    % check the algorithm options
    options = aux_checkOptions(options);

    % split the data into single data points
    points = getDataPoints(traj, true);

    % compute jacobian matrix of the dynamic function
    xsym = sym('x',[sys.nrOfDims,1]);
    usym = sym('u',[max(1,sys.nrOfInputs),1]);
    psym = sym('p',[sys.nrOfParam,1]);

    jacSym = jacobian(sys.mFile(xsym,usym,psym),psym);
    jacF = matlabFunction(jacSym,'Vars',{xsym,usym,psym});

    % optimize single-step-error using fmincon
    optOpts = optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
                            'Display','off');

    p0 = ones(sys.nrOfParam,1);

    p = fmincon(@(p) aux_costSingleStep(p,points,sys.mFile,jacF),p0, ...
                                            [],[],[],[],[],[],[],optOpts);

    % optimize multi-step-error using fmincon
    if strcmp(options.alg,'multiStep')
        
        jacSym = jacobian(sys.mFile(xsym,usym,psym),xsym);
        jacFx = matlabFunction(jacSym,'Vars',{xsym,usym,psym});

        p = fmincon(@(p) aux_costMultiStep(p,traj,sys.mFile,jacF,jacFx),p, ...
                                            [],[],[],[],[],[],[],optOpts);
    end
end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkOptions(options)
% check the algorithm settings provided by the user

    % default values
    if ~isfield(options,'alg')
        options.alg = 'singleStep';
    end

    % check user input
    if ~ismember(options.alg,{'singleStep','multiStep'})
        throw(CORAerror('CORA:wrongFieldValue','options.alg', ...
                                        {'singleStep','multiStep'}));
    end

    % check if there are any redundant options specified
    redundantOptions(options,{'alg'});
end

function [c,grad] = aux_costSingleStep(p,points,f,jacF)
% cost function for minimizing the difference between the derivative
% specified by the differential equation and the derivative numerically
% estimated from the data

    c = 0; grad = zeros(1,length(p));

    % loop over all data points
    for i = 1:size(points.x,2)

        fi = f(points.x(:,i),points.u(:,i),p);
        dfi = jacF(points.x(:,i),points.u(:,i),p);

        % compute error
        c = c + sum((points.dx(:,i) - fi).^2);

        % compute gradient of the error
        grad = grad - 2*(points.dx(:,i) - fi)' * dfi;
    end
end

function [cost,grad] = aux_costMultiStep(p,traj,f,jacFp,jacFx)
% compute the value of the cost function together with the gradient for
% minimizing the difference between the given trajectory and the prediction
% by the model.

    n = size(traj(1).x,1);

    % initialization
    cost = 0;
    grad = zeros(1,length(p));

    % loop over all trajectories
    for i = 1:length(traj)

        x = traj(i).x(:,1);
        dx = zeros(n,length(p));

        % loop over all time steps
        for j = 2:length(traj(i).t)
            
            xi = traj(i).x(:,j);
            dt = traj(i).t(j) - traj(i).t(j-1);

            if isfield(traj(i),'u')
                ui = traj(i).u(:,j-1);
            else
                ui = 0;
            end

            % compute next state using Runge-Kutta-4 
            x_prev = x;

            k1 = f(x,ui,p);
            k2 = f(x + dt/2*k1,ui,p);
            k3 = f(x + dt/2*k2,ui,p);
            k4 = f(x + dt*k3,ui,p);

            x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);

            % compute error
            cost = cost + sum((x - xi).^2);

            % compute gradient dk1/dp = dk1/dp + dk1/dp * dx/dp
            dk1 = jacFp(x_prev,ui,p) + jacFx(x_prev,ui,p)*dx;

            % compute gradient dk2/dp = dk2/dp + dk2/dx * dx/dp
            dk2 = jacFp(x_prev + dt/2*k1,ui,p) + ...
                             jacFx(x_prev + dt/2*k1,ui,p)*(dx + dt/2*dk1);

            % compute gradient dk3/dp = dk3/dp + dk3/dx * dx/dp
            dk3 = jacFp(x_prev + dt/2*k2,ui,p) + ...
                             jacFx(x_prev + dt/2*k2,ui,p)*(dx + dt/2*dk2);

            % compute gradient dk4/dp = dk4/dp + dk4/dx * dx/dp
            dk4 = jacFp(x_prev + dt*k3,ui,p) + ...
                             jacFx(x_prev + dt*k3,ui,p)*(dx + dt*dk3);

            % compute overall gradient
            dx = dx + dt/6*(dk1 + 2*dk2 + 2*dk3 + dk4);

            grad = grad + 2*(x - xi)' * dx;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
