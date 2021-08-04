classdef nonlinearSysDT < contDynamics
% nonlinearSysDT class (time-discrete nonlinear system)
%
% Syntax:  
%    obj = nonlinearSysDT(fun,dt)
%    obj = nonlinearSysDT(name,fun,dt)
%    obj = nonlinearSysDT(fun,dt,states,inputs)
%    obj = nonlinearSysDT(name,fun,dt,states,inputs)
%    obj = nonlinearSysDT(name,fun,dt,states,inputs,C)
%
% Inputs:
%    fun - function handle to the dynamic equation
%    name - name of dynamics
%    dt - sampling time
%    states - number of states
%    inputs - number of inputs
%    C - output matrix
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    f = @(x,u) [x(1) + u(1);x(2) + u(2)*cos(x(1));x(3) + u(2)*sin(x(1))];
%    dt = 0.25;
%    sys = nonlinearSysDT(f,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
%               20-March-2020 (MA, simulate random removed, now provided by inherinted class)
%               19-May-2020 (NK, changed constructor syntax)
%               02-February-2021 (MW, add switching between tensor files)
%               25-March-2021 (MA, measurement matrix added)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    mFile = [];                 % function handle dynamic equation
    jacobian = [];              % function handle jacobian matrix
    hessian = [];               % function handle hessian tensor
    thirdOrderTensor = [];      % function handle third-order tensor
    linError = [];
    dt = [];                    % sampling time
    C = [];                     % measurement matrix
end

methods
    
    % class constructor
    function obj = nonlinearSysDT(varargin)
        
         name = []; states = []; inputs = []; C = [];
        
        % parse input arguments
        if nargin == 2
            fun = varargin{1};
            dt = varargin{2};
        elseif nargin == 3
            name = varargin{1};
            fun = varargin{2};
            dt = varargin{3};
        elseif nargin == 4
            fun = varargin{1};
            dt = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
        elseif nargin == 5
            name = varargin{1};
            fun = varargin{2};
            dt = varargin{3};
            states = varargin{4};
            inputs = varargin{5};
        elseif nargin == 6
            name = varargin{1};
            fun = varargin{2};
            dt = varargin{3};
            states = varargin{4};
            inputs = varargin{5};
            C = varargin{6};
        else
            error('Wrong number of input arguments!');
        end

        % get name from function handle
        if isempty(name)    
            name = func2str(fun);
            name = replace(name,{'@','(',')',','},'');
            if ~isvarname(name)
                name = 'nonlinearSysDT';
            end
        end
        
        % get number of states and number of inputs 
        if isempty(states) || isempty(inputs)
            try
                temp = numberOfInputs(fun,2);
 
                states = temp(1);
                inputs = max(1,temp(2));
            catch
                error(['Failed to determine number of states and ' ...
                       'inputs automatically! Please provide number of ' ...
                       'states and inputs as additional input arguments!']); 
            end
        end
        
        % number of outputs (for contDynamics constructor)
        outputs = 1;
        if ~isempty(C)
        	outputs = size(C,1); 
        end
        
        % generate parent object
        obj@contDynamics(name,states,inputs,outputs);
        
        % assign object properties
        obj.mFile = fun;
        obj.dt = dt;
        obj.C = C;
        
        obj.jacobian = eval(['@jacobian_',name]);
        obj.hessian = eval(['@hessianTensor_' obj.name]);
        obj.thirdOrderTensor = eval(['@thirdOrderTensor_' obj.name]);
        
    end
    
    
    function obj = setHessian(obj,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            obj.hessian = eval(['@hessianTensor_' obj.name]);
        elseif strcmp(version,'int')
            obj.hessian = eval(['@hessianTensorInt_' obj.name]);
        end
    end
    function obj = setThirdOrderTensor(obj,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            obj.thirdOrderTensor = eval(['@thirdOrderTensor_' obj.name]);
        elseif strcmp(version,'int')
            obj.thirdOrderTensor = eval(['@thirdOrderTensorInt_' obj.name]);
        end
    end
    
end
end

%------------- END OF CODE --------------