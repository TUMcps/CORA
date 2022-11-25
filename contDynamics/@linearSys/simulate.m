function [t,x,ind] = simulate(obj,params,varargin)
% simulate - simulates a linear system 
%
% Syntax:  
%    [t,x] = simulate(obj,params)
%    [t,x,ind] = simulate(obj,params,options)
%
% Inputs:
%    obj - linearSys object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .u: input signal
%    options - ODE45 options (optional)
%
% Outputs:
%    t - time vector
%    x - state trajectory
%    ind - returns the event which has been detected
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

% Author:        Matthias Althoff
% Written:       03-May-2007 
% Last update:   20-March-2008
%                08-May-2020 (MW, update interface)
% Last revision: ---

%------------- BEGIN CODE --------------

    % parse input arguments
    isOpt = 0;
    if nargin >= 3 && ~isempty(varargin{1})
       options = varargin{1}; 
       isOpt = 1;
    end
    
    if ~isfield(params,'u')
       params.u = zeros(obj.nrOfInputs,1); 
    end
    
    if ~isfield(params,'tStart')
       params.tStart = 0;
    end
    
    tFinal = (params.tFinal-params.tStart)/size(params.u,2);
    
    if isfield(params,'timeStep')
        tSpan = 0:params.timeStep:tFinal;
        if abs(tSpan(end)-tFinal) > 1e-10
           tSpan = [tSpan,tFinal]; 
        end
    else
        tSpan = [0,tFinal];
    end

    % simulate the system
    params_ = params;
    t = [];
    x = [];
    ind = [];
    x0 = params.x0;
    
    for i = 1:size(params.u,2)
        
        params_.u = params.u(:,i);
        
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
        x = [x;x_]; 
        if isempty(t)
            t = t_ + params.tStart;
        else
            t = [t;t_ + t(end)]; 
        end
        x0 = x(end,:)';
        
        if ~isempty(ind)
           return; 
        end
    end

%------------- END OF CODE --------------