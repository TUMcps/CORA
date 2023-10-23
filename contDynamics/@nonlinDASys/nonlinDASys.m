classdef nonlinDASys < contDynamics
% nonlinDASys class (nonlinear differential algebraic system)
%    x' = f(x,y,u)    % dynamic equation
%    0  = g(x,y,u)    % constraint equation
%    z  = h(x,y,u)    % output equation
%
% Syntax:
%    % only dynamic and constraint equation
%    obj = nonlinDASys(dynFun,conFun)
%    obj = nonlinDASys(name,dynFun,conFun)
%    obj = nonlinDASys(dynFun,conFun,states,inputs,constraints)
%    obj = nonlinDASys(name,dynFun,conFun,states,inputs,constraints)
%
%    % dynamic and constraint equation with output equation
%    obj = nonlinDASys(dynFun,conFun,outFun)
%    obj = nonlinDASys(name,dynFun,conFun,outFun)
%    obj = nonlinDASys(dynFun,conFun,states,inputs,constraints,outFun,outputs)
%    obj = nonlinDASys(name,dynFun,conFun,states,inputs,constraints,outFun,outputs)
%
% Inputs:
%    name - name of the system
%    dynFun - function handle to dynamic equation
%    conFun - function handle to constraint equation
%    states - number of states
%    inputs - number of inputs
%    constraints - number of constraints
%    outFun - function handle to output equation
%    outputs - number of outputs
%
% Outputs:
%    obj - generated nonlinDASys object
%
% Example:
%    f = @(x,y,u) x(1)+1+u(1);
%    g = @(x,y,u) (x(1)+1)*y(1) + 2;
%    sys = nonlinDASys(f,g);
%
%    h = @(x,y,u) x(1) + y(1);
%    sys = nonlinDASys(f,g,h);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       27-October-2011
% Last update:   02-February-2021 (MW, add switching between tensor files)
%                18-November-2022 (MW, add output equation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % dynamic and constraint equations
    nrOfConstraints = 0;            % number of constraints
    dynFile = [];                   % function handle to dynamic file
    conFile = [];                   % function handle to constraint file
    jacobian = [];                  % function handle to jacobian matrix
    hessian = [];                   % function handle to hessian tensor
    thirdOrderTensor = [];          % function handle to third-order tensor

    % output equations
    out_mFile = [];                 % function handle to output equatio
    out_isLinear = [];              % which output functions are linear
    out_jacobian = [];              % function handle to jacobian matrix
    out_hessian = [];               % function handle to hessian tensor
    out_thirdOrderTensor = [];      % function handle to third-order tensor

    linError = [];	                % linearization error
end
    
methods
    
    % class constructor
    function obj = nonlinDASys(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'nonlinDASys')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,dynFun,conFun,states,inputs,constraints,outFun,outputs] = ...
            aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,dynFun,conFun,states,inputs,constraints,outFun,outputs);
        
        % 4. analyze functions and extract number of states, inputs, constraints, outputs
        [states,inputs,constraints,outFun,outputs,out_isLinear] = ...
            aux_computeProperties(dynFun,conFun,states,inputs,constraints,outFun,outputs);
        
        % 5. instantiate parent class and assign object properties
        obj@contDynamics(name,states,inputs,outputs);
        obj.nrOfConstraints = constraints;
        obj.dynFile = dynFun;
        obj.conFile = conFun;
        obj.out_mFile = outFun;
        obj.out_isLinear = out_isLinear;
        % link jacobian, hessian and third-order tensor files
        obj.jacobian = eval(['@jacobian_' obj.name]);
        obj.hessian = eval(['@hessianTensor_' obj.name]);
        obj.thirdOrderTensor = eval(['@thirdOrderTensor_' obj.name]);
        obj.out_jacobian = eval(['@out_jacobian_',name]);
        obj.out_hessian = eval(['@out_hessianTensor_',name]);
        obj.out_thirdOrderTensor = eval(['@out_thirdOrderTensor_',name]);
    end
    
    function obj = setHessian(obj,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            obj.hessian = eval(['@hessianTensor_' obj.name]);
        elseif strcmp(version,'int')
            obj.hessian = eval(['@hessianTensorInt_' obj.name]);
        end
    end
    function obj = setOutHessian(obj,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            obj.out_hessian = eval(['@out_hessianTensor_' obj.name]);
        elseif strcmp(version,'int')
            obj.out_hessian = eval(['@out_hessianTensorInt_' obj.name]);
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
    function obj = setOutThirdOrderTensor(obj,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            obj.out_thirdOrderTensor = eval(['@out_thirdOrderTensor_' obj.name]);
        elseif strcmp(version,'int')
            obj.out_thirdOrderTensor = eval(['@out_thirdOrderTensorInt_' obj.name]);
        end
    end
    
end
end


% Auxiliary functions -----------------------------------------------------

function [name,dynFun,conFun,states,inputs,constraints,outFun,outputs] = ...
            aux_parseInputArgs(varargin)

    % check number of input arguments
    if nargin ~= 0 && nargin < 2
        throw(CORAerror('CORA:notEnoughInputArgs',2));
    elseif nargin > 8
        throw(CORAerror('CORA:tooManyInputArgs',8));
    end

    % default values
    name = []; states = []; inputs = []; constraints = [];
    outFun = []; outputs = [];

    % no input arguments
    if nargin == 0
        return;
    end

    % parse input arguments
    if nargin == 2
        % syntax: obj = nonlinDASys(dynFun,conFun)
        dynFun = varargin{1};
        conFun = varargin{2};
    elseif nargin == 3
        if ischar(varargin{1})
            % syntax: obj = nonlinDASys(name,dynFun,conFun)
            name = varargin{1};
            dynFun = varargin{2};
            conFun = varargin{3};
        elseif isa(varargin{1},'function_handle')
            % syntax: obj = nonlinDASys(dynFun,conFun,outFun)
            dynFun = varargin{1};
            conFun = varargin{2};
            outFun = varargin{3};
        end
    elseif nargin == 4
        % syntax: obj = nonlinDASys(name,dynFun,conFun,outFun)
        name = varargin{1};
        dynFun = varargin{2};
        conFun = varargin{3};
        outFun = varargin{4};
    elseif nargin == 5
        % syntax: obj = nonlinDASys(dynFun,conFun,states,inputs,constraints)
        dynFun = varargin{1};
        conFun = varargin{2};
        states = varargin{3};
        inputs = varargin{4};
        constraints = varargin{5};
    elseif nargin == 6
        % syntax: obj = nonlinDASys(name,dynFun,conFun,states,inputs,constraints)
        name = varargin{1};
        dynFun = varargin{2};
        conFun = varargin{3};
        states = varargin{4};
        inputs = varargin{5};
        constraints = varargin{6};
    elseif nargin == 7
        % syntax: obj = nonlinDASys(dynFun,conFun,states,inputs,constraints,outFun,outputs)
        dynFun = varargin{1};
        conFun = varargin{2};
        states = varargin{3};
        inputs = varargin{4};
        constraints = varargin{5};
        outFun = varargin{6};
        outputs = varargin{7};
    elseif nargin == 8
        % syntax: obj = nonlinDASys(name,dynFun,conFun,states,inputs,constraints,outFun,outputs)
        name = varargin{1};
        dynFun = varargin{2};
        conFun = varargin{3};
        states = varargin{4};
        inputs = varargin{5};
        constraints = varargin{6};
        outFun = varargin{7};
        outputs = varargin{8};
    end
    
    % get name from function handle
    if isempty(name)    
        name = func2str(dynFun);
        name = replace(name,{'@','(',')',','},'');
        if ~isvarname(name)
            name = 'nonlinDASys';
        end
    end

end

function aux_checkInputArgs(name,dynFun,conFun,states,inputs,constraints,outFun,outputs)

    % check name (not empty because default name is not empty)
    if ~ischar(name)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'System name has to be a char array.'));
    end
    
    % fun and out_fun have to be function handles with two inputs
    if ~isempty(dynFun) && (~isa(dynFun,'function_handle') || nargin(dynFun) ~= 3)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Dynamic function has to be a function handle with three input arguments.'));
    end
    if ~isempty(conFun) && (~isa(conFun,'function_handle') || nargin(conFun) ~= 3)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Constraint function has to be a function handle with three input arguments.'));
    end
    if ~isempty(outFun) && (~isa(outFun,'function_handle') || nargin(outFun) ~= 3)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Output function has to be a function handle with three input arguments.'));
    end

    % states, inputs, and outputs have to be numeric, scalar integer > 0
    if ~isempty(states)
        inputArgsCheck({{states,'att','numeric',...
            {'positive','integer','scalar'}}});
    end
    if ~isempty(inputs)
        inputArgsCheck({{inputs,'att','numeric',...
            {'positive','integer','scalar'}}});
    end
    if ~isempty(constraints)
        inputArgsCheck({{constraints,'att','numeric',...
            {'positive','integer','scalar'}}});
    end
    if ~isempty(outputs)
        inputArgsCheck({{outputs,'att','numeric',...
            {'positive','integer','scalar'}}});
    end

end

function [states,inputs,constraints,outFun,outputs,out_isLinear] = ...
    aux_computeProperties(dynFun,conFun,states,inputs,constraints,outFun,outputs)

    % get number of states, inputs, and constraints 
    if isempty(states) || isempty(inputs) || isempty(constraints)
        try
            [temp1,states] = inputArgsLength(dynFun,3);
            temp2 = inputArgsLength(conFun,3);
            constraints = max(temp1(2),temp2(2));
            inputs = max(1,max(temp1(3),temp2(3)));
            % ensure that constraint equation does not contain more states
            % than dynamic equation
            if temp2(1) > states
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'More states in constraint equation than in dynamic equation.'));
            end
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of states and ' ...
                   'inputs automatically! Please provide number of ' ...
                   'states and inputs as additional input arguments!'])); 
        end
    end

    % default output equation and number of outputs (= states)
    if isempty(outFun)
        outFun = @(x,y,u) eye(states)*x(1:states);
        outputs = states;
        out_isLinear = true(outputs,1);
    else
        try
            [temp,out_out] = inputArgsLength(outFun,3);
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of outputs automatically!\n'...
                   'Please provide number of outputs ' ...
                   'as an additional input argument!']));
        end

        % use computed value if not provided
        if isempty(outputs)
            outputs = out_out;
        end
           
        % ensure that output equation does not contain more states than
        % dynamic equation
        if temp(1) > states
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'More states in output equation than in dynamic equation.'));
        end

        % check which output functions are linear
        out_isLinear = isFuncLinear(outFun,[states;constraints;inputs]);
    end

end

% ------------------------------ END OF CODE ------------------------------
