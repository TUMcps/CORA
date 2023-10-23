classdef nonlinParamSys < contDynamics
% nonlinParamSys class (nonlinear parametric system; parameters can be 
%    constant or vary over time)
%    x' = f(x,u,p)    % dynamic equation
%    y  = g(x,u,p)    % output equation
%
% Syntax:
%    % only dynamic equation
%    obj = nonlinParamSys(fun)
%    obj = nonlinParamSys(fun,type)
%    obj = nonlinParamSys(name,fun)
%    obj = nonlinParamSys(name,fun,type)
%    obj = nonlinParamSys(fun,states,inputs,params)
%    obj = nonlinParamSys(fun,states,inputs,params,type)
%    obj = nonlinParamSys(name,fun,states,inputs,params)
%    obj = nonlinParamSys(name,fun,states,inputs,params,type)
%
%    % dynamic equation and output equation
%    obj = nonlinParamSys(fun,out_fun)
%    obj = nonlinParamSys(fun,type,out_fun)
%    obj = nonlinParamSys(name,fun,out_fun)
%    obj = nonlinParamSys(name,fun,type,out_fun)
%    obj = nonlinParamSys(fun,states,inputs,params,out_fun,outputs)
%    obj = nonlinParamSys(fun,states,iputs,params,type,out_fun,outputs)
%    obj = nonlinParamSys(name,fun,states,inputs,params,out_fun,outputs)
%    obj = nonlinParamSys(name,fun,states,inputs,params,type,out_fun,outputs)
%
% Inputs:
%    fun - function handle to the dynamic equation
%    name - name of dynamics
%    type - time-dependency of parameters
%               'constParam' (constant parameter, default)
%               'varParam' (time-varying parameter)
%    states - number of states
%    inputs - number of inputs
%    params - number of parameters
%    out_fun - function handle to the output equation
%    outputs - number of outputs
%
% Outputs:
%    obj - generated nonlinParamSys object
%
% Example:
%    f = @(x,u,p) [x(2); p(1)*(1-x(1)^2)*x(2)-x(1)];
%    sys = nonlinParamSys('vanDerPol',f)
%
%    g = @(x,u,p) x(1) - p(1);
%    sys = nonlinParamSys('vanDerPol',f,g);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       23-September-2010
% Last update:   27-October-2011
%                16-August-2016
%                02-June-2017
%                18-May-2020 (NK, changed constructor syntax)
%                18-November-2022 (MW, add output equation)
%                23-November-2022 (MW, introduce checks, restructure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
  

properties (SetAccess = private, GetAccess = public)
    constParam = true;              % constant or time-varying parameters
    nrOfParam = 1;                  % number of parameters

    % dynamic equation
    mFile = [];                     % function handle to dynamic file
    jacobian = [];                  % function handle to jacobian matrix
    jacobian_freeParam = [];        % ?
    parametricDynamicFile = [];     % ?
    hessian = [];                   % function handle to hessian tensor
    thirdOrderTensor = [];          % function handle to third-order tensor

    % output equation
    out_mFile = [];                 % function handle output equation
    out_isLinear = [];              % which output functions are linear
    out_jacobian = [];              % function handle jacobian matrix
    out_jacobian_freeParam = [];    % ?
    out_parametricDynamicFile = []; % ?
    out_hessian = [];               % function handle hessian tensor
    out_thirdOrderTensor = [];      % function handle third-order tensor

    derivative = [];
    linError = [];
end
    
methods
    
    % class constructor
    function obj = nonlinParamSys(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'nonlinParamSys')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,fun,states,inputs,params,type,out_fun,outputs] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,fun,states,inputs,params,type,out_fun,outputs);

        % 4. analyze functions and extract number of states, inputs, params, outputs
        [states,inputs,params,out_fun,outputs,out_isLinear] = ...
            aux_computeProperties(fun,states,inputs,params,out_fun,outputs);

        % 5. instantiate parent class
        obj@contDynamics(name,states,inputs,outputs);
        
        % 6a. assign object properties: dynamic equation
        obj.nrOfParam = params;
        if strcmp(type,'varParam')
           obj.constParam = false;
        end
        obj.mFile = fun;
        obj.jacobian = eval(['@jacobian_',name]);
        obj.jacobian_freeParam = eval(['@jacobian_freeParam_',name]);
        obj.parametricDynamicFile = eval(['@parametricDynamicFile_',name]);
        obj.hessian = eval(['@hessianTensor_',name]);
        obj.thirdOrderTensor = eval(['@thirdOrderTensor_',name]);
        
        % 6b. assign object properties: output equation
        obj.out_mFile = out_fun;
        obj.out_isLinear = out_isLinear;
        obj.out_jacobian = eval(['@out_jacobian_',name]);
        obj.out_jacobian_freeParam = eval(['@out_jacobian_freeParam_',name]);
        obj.out_parametricDynamicFile = eval(['@out_parametricDynamicFile_',name]);
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

function [name,fun,states,inputs,params,type,out_fun,outputs] = aux_parseInputArgs(varargin)

    % check number of input arguments
    if nargin > 8
        throw(CORAerror('CORA:tooManyInputArgs',8));
    end

    % default values
    name = []; states = []; inputs = []; params = []; 
    type = 'constParam';
    out_fun = []; outputs = [];

    % no input arguments
    if nargin == 0
        return;
    end

    % parse input arguments (order for each number not fixed)
    if nargin == 1
        % syntax: obj = nonlinParamSys(fun)
        fun = varargin{1};
    elseif nargin == 2
        if ischar(varargin{1})
            % syntax: obj = nonlinParamSys(name,fun)
            name = varargin{1};
            fun = varargin{2};
        elseif ischar(varargin{2})
            % syntax: obj = nonlinParamSys(fun,type)
            fun = varargin{1};
            type = varargin{2};
        elseif isa(varargin{2},'function_handle')
            % syntax: obj = nonlinParamSys(fun,out_fun)
            fun = varargin{1};
            out_fun = varargin{2};
        end
    elseif nargin == 3
        if ischar(varargin{1})
            if ischar(varargin{3})
                % syntax: obj = nonlinParamSys(name,fun,type)
                name = varargin{1};
                fun = varargin{2};
                type = varargin{3};
            else
                % syntax: obj = nonlinParamSys(name,fun,out_fun)
                name = varargin{1};
                fun = varargin{2};
                out_fun = varargin{3};
            end
        else
            % syntax: obj = nonlinParamSys(fun,type,out_fun)
            fun = varargin{1};
            type = varargin{2};
            out_fun = varargin{3};
        end
    elseif nargin == 4
        if ischar(varargin{1})
            % syntax: obj = nonlinParamSys(name,fun,type,out_fun)
            name = varargin{1};
            fun = varargin{2};
            type = varargin{3};
            out_fun = varargin{4};
        else
            % syntax: obj = nonlinParamSys(fun,states,inputs,params)
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
            params = varargin{4};
        end
    elseif nargin == 5
        if ischar(varargin{1})
            % syntax: obj = nonlinParamSys(name,fun,states,inputs,params)
            name = varargin{1};
            fun = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
            params = varargin{5};
        else
            % syntax: obj = nonlinParamSys(fun,states,inputs,params,type)
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
            params = varargin{4};
            type = varargin{5};
        end
    elseif nargin == 6
        if ischar(varargin{1})
            % syntax: obj = nonlinParamSys(name,fun,states,inputs,params,type)
            name = varargin{1};
            fun = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
            params = varargin{5};
            type = varargin{6};
        else
            % syntax: obj = nonlinParamSys(fun,states,inputs,params,out_fun,outputs)
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
            params = varargin{4};
            out_fun = varargin{5};
            outputs = varargin{6};
        end
    elseif nargin == 7
        if ischar(varargin{1})
            % syntax: obj = nonlinParamSys(name,fun,states,inputs,params,out_fun,outputs)
            name = varargin{1};
            fun = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
            params = varargin{5};
            out_fun = varargin{6};
            outputs = varargin{7};
        else
            % syntax: obj = nonlinParamSys(fun,states,inputs,params,type,out_fun,outputs)
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
            params = varargin{4};
            type = varargin{5};
            out_fun = varargin{6};
            outputs = varargin{7};
        end
    elseif nargin == 8
        % syntax: obj = nonlinParamSys(name,fun,states,inputs,params,type,out_fun,outputs)
        name = varargin{1};
        fun = varargin{2};
        states = varargin{3};
        inputs = varargin{4};
        params = varargin{5};
        type = varargin{6};
        out_fun = varargin{7};
        outputs = varargin{8};
    end

    % get name from function handle
    if isempty(name)    
        name = func2str(fun);
        name = replace(name,{'@','(',')',','},'');
        if ~isvarname(name)
            name = 'nonlinParamSys';
        end
    end

end

function aux_checkInputArgs(name,fun,states,inputs,params,type,out_fun,outputs)

    if CHECKS_ENABLED
    % check name (not empty because default name is not empty)
    if ~ischar(name)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'System name has to be a char array.'));
    end
    
    % fun and out_fun have to be function handles with two inputs
    if ~isempty(fun) && (~isa(fun,'function_handle') || nargin(fun) ~= 3)
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Dynamic function has to be a function handle with three input arguments.'));
    end
    if ~isempty(out_fun) && (~isa(out_fun,'function_handle') || nargin(out_fun) ~= 3)
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
    if ~isempty(params)
        inputArgsCheck({{params,'att','numeric',...
            {'positive','integer','scalar'}}});
    end
    if ~isempty(outputs)
        inputArgsCheck({{outputs,'att','numeric',...
            {'positive','integer','scalar'}}});
    end
    
    % char has a default value (so cannot be empty)
    if ~ischar(type) || ~ismember(type,{'constParam','varParam'})
        throw(CORAerror('CORA:wrongInputInConstructor',...
            '"Type" has to be either "constParam" or "varParam".'));
    end
    end

end

function [states,inputs,params,out_fun,outputs,out_isLinear] = ...
            aux_computeProperties(fun,states,inputs,params,out_fun,outputs)

    % get number of states and number of inputs 
    if isempty(states) || isempty(inputs) || isempty(params)
        try
            [temp,states] = inputArgsLength(fun,3);
            % CORA models have to have at least one input
            inputs = max(1,temp(2));
            params = max(1,temp(3));
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of states and ' ...
                   'inputs automatically! Please provide number of ' ...
                   'states, inputs, and parameter as additional ', ...
                   'input arguments!'])); 
        end
    end

    % default output equation and number of outputs (= states)
    if isempty(out_fun)
        out_fun = @(x,u,p) eye(states)*x(1:states);
        outputs = states;
        out_isLinear = true(outputs,1);

    else
        try
            [temp,out_out] = inputArgsLength(out_fun,3);
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
        out_isLinear = isFuncLinear(out_fun,[states;inputs;outputs]);
    end
end

% ------------------------------ END OF CODE ------------------------------
