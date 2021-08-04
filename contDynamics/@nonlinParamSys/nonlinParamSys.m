classdef nonlinParamSys < contDynamics
% nonlinParamSys class (nonlinear parametric system; parameters can be 
% constant or vary over time)
%
% Syntax:  
%    obj = nonlinParamSys(fun)
%    obj = nonlinParamSys(fun,type)
%    obj = nonlinParamSys(name,fun)
%    obj = nonlinParamSys(name,fun,type)
%    obj = nonlinParamSys(fun,states,inputs,params)
%    obj = nonlinParamSys(fun,states,iputs,params,type)
%    obj = nonlinParamSys(name,fun,states,inputs,params)
%    obj = nonlinParamSys(name,fun,states,inputs,params,type)
%
% Inputs:
%    fun - function handle to the dynamic equation
%    name - name of dynamics
%    type - 'constParam' (constant parameter, default) or 'varParam' (time
%            varying parameter)
%    states - number of states
%    inputs - number of inputs
%    params - number of parameter
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    f = @(x,u,p) [x(2); ...
%                  p(1)*(1-x(1)^2)*x(2)-x(1)];
%
%    sys = nonlinParamSys('vanDerPol',f)
%
% Outputs:
%    obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      23-September-2010
% Last update:  27-October-2011
%               16-August-2016
%               02-June-2017
%               18-May-2020 (NK, changed constructor syntax)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    nrOfParam = 1;
    mFile = [];
    jacobian = [];
    jacobian_freeParam = [];
    hessian = [];
    thirdOrderTensor = [];
    parametricDynamicFile = [];
    derivative = [];
    linError = [];
    constParam = 1;
end
    
methods
    
    % class constructor
    function obj = nonlinParamSys(varargin)
        
        name = []; states = []; inputs = []; params = []; 
        type = 'constParam';
        
        % parse input arguments
        if nargin == 1
            fun = varargin{1};
        elseif nargin == 2
            if ischar(varargin{1})
                name = varargin{1};
                fun = varargin{2};
            else
                fun = varargin{1};
                type = varargin{2};
            end
        elseif nargin == 3
            name = varargin{1};
            fun = varargin{2};
            type = varargin{3};
        elseif nargin == 4
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
            params = varargin{4};
        elseif nargin == 5
            if ischar(varargin{1})
                name = varargin{1};
                fun = varargin{2};
                states = varargin{3};
                inputs = varargin{4};
                params = varargin{5};
            else
                fun = varargin{1};
                states = varargin{2};
                inputs = varargin{3};
                params = varargin{4};
                type = varargin{5};
            end
        elseif nargin == 6
                name = varargin{1};
                fun = varargin{2};
                states = varargin{3};
                inputs = varargin{4};
                params = varargin{5};
                type = varargin{6};
        else
            error('Wrong number of input arguments!');
        end

        % check input arguments
        if ~ischar(type) || ~ismember(type,{'constParam','varParam'})
           error('Wrong value for input argument "type"!'); 
        end
        
        % get name from function handle
        if isempty(name)    
            name = func2str(fun);
            name = replace(name,{'@','(',')',','},'');
            if ~isvarname(name)
                name = 'nonlinParamSys';
            end
        end
        
        % get number of states and number of inputs 
        if isempty(states) || isempty(inputs) || isempty(params)
            try
                [temp,states] = numberOfInputs(fun,3);
                inputs = max(1,temp(2));
                params = max(1,temp(3));
            catch
                error(['Failed to determine number of states and ' ...
                       'inputs automatically! Please provide number of ' ...
                       'states, inputs, and parameter as additional ', ...
                       'input arguments!']); 
            end
        end
        
        % generate parent object
        obj@contDynamics(name,states,inputs,1);
        
        % assign object properties
        obj.nrOfParam = params;
        obj.mFile = fun;
        
        if strcmp(type,'varParam')
           obj.constParam = 0; 
        end

        obj.jacobian = eval(['@jacobian_',name]);
        obj.jacobian_freeParam = eval(['@jacobian_freeParam_',name]);
        obj.hessian = eval(['@hessianTensor_',name]);
        obj.thirdOrderTensor = eval(['@thirdOrderTensor_',name]);
        obj.parametricDynamicFile = eval(['@parametricDynamicFile_',name]);
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