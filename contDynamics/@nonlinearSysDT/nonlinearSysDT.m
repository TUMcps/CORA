classdef nonlinearSysDT < contDynamics
% nonlinearSysDT class (time-discrete nonlinear system)
%    x_k+1 = f(x_k,u_k)
%    y_k   = g(x_k,u_k)
%
% Syntax:
%    % only dynamic equation
%    obj = nonlinearSysDT(fun,dt)
%    obj = nonlinearSysDT(name,fun,dt)
%    obj = nonlinearSysDT(fun,dt,states,inputs)
%    obj = nonlinearSysDT(name,fun,dt,states,inputs)
%
%    % dynamic equation and output equation
%    obj = nonlinearSysDT(fun,dt,out_fun)
%    obj = nonlinearSysDT(name,fun,dt,out_fun)
%    obj = nonlinearSysDT(fun,dt,states,inputs,out_fun,outputs)
%    obj = nonlinearSysDT(name,fun,dt,states,inputs,out_fun,outputs)
%
% Inputs:
%    fun - function handle to the dynamic equation
%    name - name of dynamics
%    dt - sampling time
%    states - number of states
%    inputs - number of inputs
%    out_fun - function handle to the output equation
%    outputs - number of outputs
%
% Outputs:
%    obj - generated nonlinearSysDT object
%
% Example:
%    f = @(x,u) [x(1) + u(1);x(2) + u(2)*cos(x(1));x(3) + u(2)*sin(x(1))];
%    dt = 0.25;
%    sys = nonlinearSysDT(f,dt)
%
%    g = @(x,u) x(1) + x(2);
%    sys = nonlinearSysDT(f,dt,g)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       21-August-2012
% Last update:   29-January-2018
%                20-March-2020 (MA, simulate random removed, now provided by inherinted class)
%                19-May-2020 (NK, changed constructor syntax)
%                02-February-2021 (MW, add switching between tensor files)
%                25-March-2021 (MA, measurement matrix added)
%                18-November-2022 (MW, add output equation)
%                26-June-2023 (LL, support of 2D out_isLinear-array)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
  

properties (SetAccess = private, GetAccess = public)
    % dynamic equation
    mFile = [];                 % function handle dynamic equation
    jacobian = [];              % function handle jacobian matrix
    hessian = [];               % function handle hessian tensor
    thirdOrderTensor = [];      % function handle third-order tensor
    dt {mustBeNumeric} = [];    % sampling time
    
    % output equation
    C = [];                     % measurement matrix (used to rewrite
                                % output function handle as linear equation)
    out_mFile = [];             % function handle output equation
    out_isLinear = [];          % which output functions are linear
    out_jacobian = [];          % function handle jacobian matrix
    out_hessian = [];           % function handle hessian tensor
    out_thirdOrderTensor = [];  % function handle third-order tensor

    linError = [];              % linearization error
end

methods
    
    % class constructor
    function obj = nonlinearSysDT(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
        %         if nargin == 1 && isa(varargin{1},'nonlinearSysDT')
        %             obj = varargin{1}; return
        %         end

        % 2. parse input arguments: varargin -> vars
        [name,fun,dt,states,inputs,out_fun,outputs] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,fun,dt,states,inputs,out_fun,outputs);

        % 4. default output equation and number of outputs (= states)
        [states,inputs,out_fun,outputs,out_isLinear,rewriteAsC,C] = ...
            aux_computeProperties(fun,states,inputs,out_fun,outputs);
        
        % 5. instantiate parent class
        obj@contDynamics(name,states,inputs,outputs);
        
        % 6a. assign object properties: dynamic equation
        obj.dt = dt;
        
        obj.mFile = fun;
        obj.jacobian = eval(['@jacobian_',name]);
        obj.hessian = eval(['@hessianTensor_' obj.name]);
        obj.thirdOrderTensor = eval(['@thirdOrderTensor_' obj.name]);

        % 6b. assign object properties: output equation
        obj.out_mFile = out_fun;
        obj.out_isLinear = out_isLinear;
        if all(all(out_isLinear)) && rewriteAsC
            C = aux_rewriteOutFunAsMatrix(out_fun,states,outputs);
        end
        obj.C = C;
        obj.out_jacobian = eval(['@out_jacobian_',name]);
        obj.out_hessian = eval(['@out_hessianTensor_',name]);
        obj.out_thirdOrderTensor = eval(['@out_thirdOrderTensor_',name]);
        
    end
    
    % set tensors to either numeric or interval arithmetic
    % (required in computation of Lagrange remainder)
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

function [name,fun,dt,states,inputs,out_fun,outputs] = aux_parseInputArgs(varargin)

    if nargin ~= 0 && nargin < 2
        throw(CORAerror('CORA:notEnoughInputArgs',2));
    elseif nargin > 7
        throw(CORAerror('CORA:tooManyInputArgs',7));
    end

    % default values
    name = []; states = []; inputs = [];
    out_fun = []; outputs = [];

    % no input arguments
    if nargin == 0
        return;
    end

    % parse input arguments
    if nargin == 2
        % syntax: obj = nonlinearSysDT(fun,dt)
        fun = varargin{1};
        dt = varargin{2};
    elseif nargin == 3
        if ischar(varargin{1})
            % syntax: obj = nonlinearSysDT(name,fun,dt)
            name = varargin{1};
            fun = varargin{2};
            dt = varargin{3};
        elseif isa(varargin{1},'function_handle')
	        % syntax: obj = nonlinearSysDT(fun,dt,out_fun)
            fun = varargin{1};
            dt = varargin{2};
            out_fun = varargin{3};
        end
    elseif nargin == 4
        if ischar(varargin{1})
            % syntax: obj = nonlinearSysDT(name,fun,dt,out_fun)
            name = varargin{1};
            fun = varargin{2};
            dt = varargin{3};
            out_fun = varargin{4};
        elseif isa(varargin{1},'function_handle')
            % syntax: obj = nonlinearSysDT(fun,dt,states,inputs)
            fun = varargin{1};
            dt = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
        end
    elseif nargin == 5
        % syntax: obj = nonlinearSysDT(name,fun,dt,states,inputs)
        name = varargin{1};
        fun = varargin{2};
        dt = varargin{3};
        states = varargin{4};
        inputs = varargin{5};
    elseif nargin == 6
        % syntax: obj = nonlinearSysDT(fun,dt,states,inputs,out_fun,outputs)
        fun = varargin{1};
        dt = varargin{2};
        states = varargin{3};
        inputs = varargin{4};
        out_fun = varargin{5};
        outputs = varargin{6};
    elseif nargin == 7
        % syntax: obj = nonlinearSysDT(name,fun,dt,states,inputs,out_fun,outputs)
        name = varargin{1};
        fun = varargin{2};
        dt = varargin{3};
        states = varargin{4};
        inputs = varargin{5};
        out_fun = varargin{6};
        outputs = varargin{7};
    end

    % get name from function handle
    if isempty(name)    
        name = func2str(fun);
        name = replace(name,{'@','(',')',','},'');
        if ~isvarname(name)
            name = 'nonlinearSysDT';
        end
    end
end

function aux_checkInputArgs(name,fun,dt,states,inputs,out_fun,outputs)

    % check name (not empty because default name is not empty)
    if ~ischar(name)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'System name has to be a char array.'));
    end

    % sampling time has to be a scalar larger than zero
    if ~isempty(dt)
        inputArgsCheck({{dt,'att','numeric',{'positive','scalar'}}});
    end
    
    % fun and out_fun have to be function handles with two inputs
    if ~isempty(fun) && (~isa(fun,'function_handle') || nargin(fun) ~= 2)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Dynamic function has to be a function handle with two input arguments.'));
    end
    if ~isempty(out_fun) && (~isa(out_fun,'function_handle') || nargin(out_fun) ~= 2)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Output function has to be a function handle with two input arguments.'));
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
    if ~isempty(outputs)
        inputArgsCheck({{outputs,'att','numeric',...
            {'positive','integer','scalar'}}});
    end

end

function C = aux_rewriteOutFunAsMatrix(out_fun,states,outputs)      

    % initialize symbolic state variables
    x = sym('x',[states,1]);

    % initaliaze output matrix
    C = zeros(outputs,states);

    % loop over each column (=number of states)
    for j=1:states

        % set other symbolic variables to 0
        x_temp = subs(x,x(1:j-1),zeros(j-1,1));
        x_temp = subs(x_temp,x(j+1:end),zeros(states-j,1));

        % evaluate output equation
        out_lhs = out_fun(x_temp);

        % number of each row (=number of outputs)
        for i=1:outputs
            if has(out_lhs(i),x(j))
                % divide and extract remaining constant
                C(i,j) = double(out_lhs(i)/x(j));
            end
        end
    end

end

function [states,inputs,out_fun,outputs,out_isLinear,rewriteAsC,C] = ...
        aux_computeProperties(fun,states,inputs,out_fun,outputs)

    % get number of states and number of inputs 
    if isempty(states) || isempty(inputs)
        try
            temp = inputArgsLength(fun,2);
            states = temp(1);
            inputs = max(1,temp(2));
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of states and ' ...
                   'inputs automatically! Please provide number of ' ...
                   'states and inputs as additional input arguments!'])); 
        end
    end

    % init linear output matrix
    C = [];

    if isempty(out_fun)
        out_fun = @(x,u) eye(states)*x(1:states);
        outputs = states;

        % all output equations are linear
        out_isLinear = true(outputs,1);

        % define equivalent linear output matrix
        C = eye(states);

        % set flag to false to avoid recomputation
        rewriteAsC = false;

    else
        % compute number of outputs and number of inputs to the output
        % function (required for potential rewriting into obj.C)
        % number of inputs to the output functions
        try
            [temp,out_out] = inputArgsLength(out_fun,2);
            % number of inputs to output function has to be zero in order
            % to potentially rewrite output equation to a C matrix
            out_inputs = temp(2);
            rewriteAsC = out_inputs == 0;
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of outputs automatically!\n'...
                   'Please provide number of outputs ' ...
                   'as an additional input argument!']));
        end

        % ensure that output function does not use too many states
        if temp(1) > states
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'More states in output equation than in dynamic equation.'));
        end

        % take computed value if not provided
        if isempty(outputs)
            outputs = out_out;
        end

        % check which output functions are linear
        out_isLinear = isFuncLinear(out_fun,[states;out_inputs]);
    end

end

% ------------------------------ END OF CODE ------------------------------
