classdef nonlinDASys < contDynamics
% nonlinDASys class (nonlinear differential algebraic system)
%
% Syntax:  
%    obj = nonlinDASys(dynFun,conFun)
%    obj = nonlinDASys(name,dynFun,conFun)
%    obj = nonlinDASys(dynFun,conFun,states,inputs,constraints)
%    obj = nonlinDASys(name,dynFun,conFun,states,inputs,constraints)
%
% Inputs:
%    dynFun - function handle to dynamic equation
%    conFun - function handle to constraint equation
%    name - name of the system
%    states - number of states
%    inputs - number of inputs
%    constraints - number of constraints
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    f = @(x,y,u) x + 1 + u;
%    g = @(x,y,u) (x+1)*y + 2;
%
%    sys = nonlinDASys(f,g)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      27-October-2011
% Last update:  02-February-2021 (MW, add switching between tensor files)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    nrOfConstraints = 0;        % number of constraints
    dynFile = [];               % function handle to dynamic file
    conFile = [];               % function handle to constraint file
    jacobian = [];              % function handle to jacobian matrix
    hessian = [];               % function handle to hessian tensor
    thirdOrderTensor = [];      % function handle to third-order tensor
    linError = [];
end
    
methods
    
    % class constructor
    function obj = nonlinDASys(varargin)
        
        name = []; states = []; inputs = []; constraints = [];
        
        % parse input arguments
        if nargin == 2
            dynFun = varargin{1};
            conFun = varargin{2};
        elseif nargin == 3
            name = varargin{1};
            dynFun = varargin{2};
            conFun = varargin{3};
        elseif nargin == 5
            dynFun = varargin{1};
            conFun = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
            constraints = varargin{5};
        elseif nargin == 6
            name = varargin{1};
            dynFun = varargin{2};
            conFun = varargin{3};
            states = varargin{4};
            inputs = varargin{5};
            constraints = varargin{6};
        else
            error('Wrong number of input arguments!');
        end
        
        % get name from function handle
        if isempty(name)    
            name = func2str(dynFun);
            name = replace(name,{'@','(',')',','},'');
            if ~isvarname(name)
                name = 'nonlinDASys';
            end
        end
        
        % get number of states, inputs, and constraints 
        if isempty(states) || isempty(inputs) || isempty(constraints)
            try
                [temp1,states] = numberOfInputs(dynFun,3);
                temp2 = numberOfInputs(conFun,3);
                constraints = max(temp1(2),temp2(2));
                inputs = max(1,max(temp1(3),temp2(3)));
            catch
                error(['Failed to determine number of states and ' ...
                       'inputs automatically! Please provide number of ' ...
                       'states and inputs as additional input arguments!']); 
            end
        end
        
        % generate parent object
        obj@contDynamics(name,states,inputs,1);
        
        % assign object properties
        obj.nrOfConstraints = constraints;
        obj.dynFile = dynFun;
        obj.conFile = conFun;

        % link jacobian, hessian and third-order tensor files
        obj.jacobian = eval(['@jacobian_' obj.name]);
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