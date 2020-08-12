classdef nonlinearSys < contDynamics
% nonlinearSys class (time-continuous nonlinear system)
%
% Syntax:  
%    obj = nonlinearSys(fun)
%    obj = nonlinearSys(name,fun)
%    obj = nonlinearSys(fun,states,inputs)
%    obj = nonlinearSys(name,fun,states,inputs)
%
% Inputs:
%    fun - function handle to the dynamic equation
%    name - name of dynamics
%    states - number of states
%    inputs - number of inputs
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    f = @(x,u) [x(2); ...
%               (1-x(1)^2)*x(2)-x(1)];
%
%    sys = nonlinearSys('vanDerPol',f)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      17-October-2007 
% Last update:  29-October-2007
%               04-August-2016 (changed to new OO format)
%               19-May-2020 (NK, changed constructor syntax)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    mFile = [];                 % function handle dynamic equation
    jacobian = [];              % function handle jacobian matrix
    hessian = [];               % function handle hessian tensor
    thirdOrderTensor = [];      % function handle third-order tensor
    tensors = [];               % function handle higher order tensors
    linError = [];
end

methods
    
    % class constructor
    function obj = nonlinearSys(varargin)
        
        name = []; states = []; inputs = [];
        
        % parse input arguments
        if nargin == 1
            fun = varargin{1};
        elseif nargin == 2
            name = varargin{1};
            fun = varargin{2};
        elseif nargin == 3
            fun = varargin{1};
            states = varargin{2};
            inputs = varargin{3};
        elseif nargin == 4
            name = varargin{1};
            fun = varargin{2};
            states = varargin{3};
            inputs = varargin{4};
        else
            error('Wrong number of input arguments!');
        end

        % get name from function handle
        if isempty(name)    
            name = func2str(fun);
            name = strrep(name,'@',''); 
            name = strrep(name,'(',''); 
            name = strrep(name,')',''); 
            name = strrep(name,',','');
            if ~isvarname(name)
                name = 'nonlinearSys';
            end
        end
        
        % get number of states and number of inputs 
        if isempty(states) || isempty(inputs)
            try
                [temp,states] = numberOfInputs(fun,2);
                inputs = max(1,temp(2));
            catch
                error(['Failed to determine number of states and ' ...
                       'inputs automatically! Please provide number of ' ...
                       'states and inputs as additional input arguments!']); 
            end
        end
        
        % generate parent object
        obj@contDynamics(name,states,inputs,1);
        
        % assign object properties
        obj.mFile = fun;
        
        str = ['obj.jacobian = @jacobian_',name,';'];
        eval(str);
        
        str = ['obj.hessian = @hessianTensor_',name,';'];
        eval(str);
        
        str = ['obj.thirdOrderTensor = @thirdOrderTensor_',name,';'];
        eval(str);
        
        for i = 4:10
            str = sprintf('obj.tensors{%i} = @tensor%i_%s;',i-3,i,name);
            eval(str);
        end
    end
    
end
end

%------------- END OF CODE --------------