classdef nonlinearSys < contDynamics
% nonlinearSys class (continuous-time nonlinear system):
%    x' = f(x,u)    % dynamic equation
%    y  = g(x,u)    % output equation
%
% Syntax:
%    % only dynamic equation
%    nlnsys = nonlinearSys(fun)
%    nlnsys = nonlinearSys(name,fun)
%    nlnsys = nonlinearSys(fun,states,inputs)
%    nlnsys = nonlinearSys(name,fun,states,inputs)
%
%    % dynamic equation and output equation
%    nlnsys = nonlinearSys(fun,out_fun)
%    nlnsys = nonlinearSys(name,fun,out_fun)
%    nlnsys = nonlinearSys(fun,states,inputs,out_fun,outputs)
%    nlnsys = nonlinearSys(name,fun,states,inputs,out_fun,outputs)
%
% Inputs:
%    name - name of system
%    fun - function handle to the dynamic equation
%    states - number of states
%    inputs - number of inputs
%    out_fun - function handle to the output equation
%    outputs - number of outputs
%
% Outputs:
%    nlnsys - generated nonlinearSys object
%
% Example:
%    fun = @(x,u) [x(2); ...
%               (1-x(1)^2)*x(2)-x(1)];
%    sys = nonlinearSys('vanDerPol',fun)
%
%    out_fun = @(x,u) [x(1) + x(2)];
%    sys = nonlinearSys('vanDerPol',fun,out_fun);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       17-October-2007 
% Last update:   29-October-2007
%                04-August-2016 (changed to new OO format)
%                19-May-2020 (NK, changed constructor syntax)
%                02-February-2021 (MW, add switching between tensor files)
%                17-November-2022 (MW, add output equation)
%                23-November-2022 (MW, introduce checks, restructure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % dynamic equation
    mFile = [];                     % function handle dynamic equation
    jacobian = [];                  % function handle jacobian matrix
    hessian = [];                   % function handle hessian tensor
    thirdOrderTensor = [];          % function handle third-order tensor
    tensors = [];                   % function handle higher-order tensors

    % output equation
    out_mFile = [];                 % function handle output equation
    out_isLinear = [];              % which output functions are linear
    out_jacobian = [];              % function handle jacobian matrix
    out_hessian = [];               % function handle hessian tensor
    out_thirdOrderTensor = [];      % function handle third-order tensor

    linError = [];                  % linearization error
end

methods
    
    % class constructor
    function nlnsys = nonlinearSys(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:6,nargin);

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'nonlinearSys')
%             obj = varargin{1}; return
%         end
        
        % 2. parse input arguments: varargin -> vars
        [name,fun,states,inputs,out_fun,outputs] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,fun,states,inputs,out_fun,outputs);

        % 4. analyze functions and extract number of states, inputs, outputs
        [states,inputs,out_fun,outputs,out_isLinear] = ...
            aux_computeProperties(fun,states,inputs,out_fun,outputs);
        
        % 5. instantiate parent class
        % note: currently, we only support unit disturbance matrices
        %       (same as number of states) and unit noise matrices (same as
        %       number of outputs)
        nlnsys@contDynamics(name,states,inputs,outputs,states,outputs);
        
        % 6a. assign object properties: dynamic equation
        nlnsys.mFile = fun;
        nlnsys.jacobian = eval(['@jacobian_',name]);
        nlnsys.hessian = eval(['@hessianTensor_',name]);
        nlnsys.thirdOrderTensor = eval(['@thirdOrderTensor_',name]);
        for i = 4:10
            nlnsys.tensors{i-3} = eval(sprintf('@tensor%i_%s;',i,name));
        end

        % 6b. assign object properties: output equation
        nlnsys.out_mFile = out_fun;
        nlnsys.out_isLinear = out_isLinear;
        nlnsys.out_jacobian = eval(['@out_jacobian_',name]);
        nlnsys.out_hessian = eval(['@out_hessianTensor_',name]);
        nlnsys.out_thirdOrderTensor = eval(['@out_thirdOrderTensor_',name]);
        
    end
    
    function nlnsys = setHessian(nlnsys,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            nlnsys.hessian = eval(['@hessianTensor_' nlnsys.name]);
        elseif strcmp(version,'int')
            nlnsys.hessian = eval(['@hessianTensorInt_' nlnsys.name]);
        end
    end
    function nlnsys = setOutHessian(nlnsys,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            nlnsys.out_hessian = eval(['@out_hessianTensor_' nlnsys.name]);
        elseif strcmp(version,'int')
            nlnsys.out_hessian = eval(['@out_hessianTensorInt_' nlnsys.name]);
        end
    end

    function nlnsys = setThirdOrderTensor(nlnsys,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            nlnsys.thirdOrderTensor = eval(['@thirdOrderTensor_' nlnsys.name]);
        elseif strcmp(version,'int')
            nlnsys.thirdOrderTensor = eval(['@thirdOrderTensorInt_' nlnsys.name]);
        end
    end
    function nlnsys = setOutThirdOrderTensor(nlnsys,version)
        % allow switching between standard and interval arithmetic
        if strcmp(version,'standard')
            nlnsys.out_thirdOrderTensor = eval(['@out_thirdOrderTensor_' nlnsys.name]);
        elseif strcmp(version,'int')
            nlnsys.out_thirdOrderTensor = eval(['@out_thirdOrderTensorInt_' nlnsys.name]);
        end
    end
    
end

methods (Static = true)
    linsys = identify(varargin)       % identifies linear system from data
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,fun,states,inputs,out_fun,outputs] = aux_parseInputArgs(varargin)

    % default values
    name = []; fun =  @(x,u)[]; states = []; inputs = [];
    out_fun = []; outputs = [];

    % no input arguments
    if nargin == 0
        return;
    end

    % parse input arguments
    if nargin == 1
        % syntax: obj = nonlinearSys(fun)
        fun = varargin{1};
    elseif nargin == 2
        if ischar(varargin{1})
            % syntax: obj = nonlinearSys(name,fun)
            [name,fun] = varargin{:};
        elseif isa(varargin{1},'function_handle')
            % syntax: obj = nonlinearSys(fun,out_fun)
            [fun,out_fun] = varargin{:};
        end
    elseif nargin == 3
        if ischar(varargin{1})
            % syntax: obj = nonlinearSys(name,fun,out_fun)
            [name,fun,out_fun] = varargin{:};
        elseif isa(varargin{1},'function_handle')
            % syntax: obj = nonlinearSys(fun,states,inputs)
            [fun,states,inputs] = varargin{:};
        end
    elseif nargin == 4
        % syntax: obj = nonlinearSys(name,fun,states,inputs)
        [name,fun,states,inputs] = varargin{:};
    elseif nargin == 5
        % syntax: obj = nonlinearSys(fun,states,inputs,out_fun,outputs)
        [fun,states,inputs,out_fun,outputs] = varargin{:};
    elseif nargin == 6
        % syntax: obj = nonlinearSys(name,fun,states,inputs,out_fun,outputs)
        [name,fun,states,inputs,out_fun,outputs] = varargin{:};
    end

    % get name from function handle
    if isempty(name)
        name = func2str(fun);
        name = replace(name,{'@','(',')',','},'');
        if ~isvarname(name)
            % default name
            name = 'nonlinearSys';
        end
    end

end

function aux_checkInputArgs(name,fun,states,inputs,out_fun,outputs)
    
    if CHECKS_ENABLED
        % check name (not empty because default name is not empty)
        if ~ischar(name)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'System name has to be a char array.'));
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
    
        % states and outputs have to be numeric, scalar integer > 0,
        % inputs can be 0 (e.g., in parallel hybrid automata with only local
        % outputs = inputs and no global inputs)
        if ~isempty(states)
            inputArgsCheck({{states,'att','numeric',...
                {'positive','integer','scalar'}}});
        end
        if ~isempty(inputs)
            inputArgsCheck({{inputs,'att','numeric',...
                {'nonnegative','integer','scalar'}}});
        end
        if ~isempty(outputs)
            inputArgsCheck({{outputs,'att','numeric',...
                {'positive','integer','scalar'}}});
        end
    end

end

function [states,inputs,out_fun,outputs,out_isLinear] = ...
    aux_computeProperties(fun,states,inputs,out_fun,outputs)

    if ~isempty(inputs) && inputs == 0
        % CORA requires at least one input so that the internal
        % computations execute properly; some system do not have an input,
        % so it would be correct to explicitly state '0'
        inputs = 1;
    end

    % get number of states and number of inputs
    if isempty(states) || isempty(inputs)
        try
            [temp,states] = inputArgsLength(fun,2);
            % CORA models have to have at least one input
            inputs = max(1,temp(2));
        catch
            throw(CORAerror('CORA:specialError',...
                ['Failed to determine number of states and ' ...
                   'inputs automatically!\n'...
                   'Please provide number of states and inputs ' ...
                   'as additional input arguments!']));
        end
    end

    if isempty(out_fun)
        % init out_fun via eval to have numeric values 
        % within function handle
        outputs = states;
        out_fun = eval(sprintf('@(x,u) eye(%i)*x(1:%i)',outputs,outputs));
        out_isLinear = true(outputs,1);
    else
        % get number of states in output equation and outputs
        try
            [temp,out_out] = inputArgsLength(out_fun,2);
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

        % ensure that output equation does not have more states than
        % dynamic equation
        if temp(1) > states
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'More states in output equation than in dynamic equation.'));
        end
        
        % check which output functions are linear
        out_isLinear = isFuncLinear(out_fun,[states;inputs]);
    end
end

% ------------------------------ END OF CODE ------------------------------
