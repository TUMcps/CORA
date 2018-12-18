classdef nonlinearSys < contDynamics
% nonlinearSys class (nonlinearSys: nonlinear system)
%
% Syntax:  
%    object constructor: Obj = nonlinearSys(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    name - name of dynamics
%    dimension - system dimension
%    nrOfInputs - number of inputs
%    mFile - mfile of the state space equations
%    options - options struct
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      17-October-2007 
% Last update:  29-October-2007
%               04-August-2016 (changed to new OO format)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    mFile = [];
    jacobian = [];
    hessian = [];
    thirdOrderTensor = [];
    tensors = [];
    linError = [];
end

methods
    %class constructor
    function obj = nonlinearSys(varargin)
        % obtain function name
        if nargin>=3
            func_name = func2str(varargin{3});
            func_name = strrep(func_name,'@',''); %remove @
            func_name = strrep(func_name,'(',''); %remove (
            func_name = strrep(func_name,')',''); %remove )
            func_name = strrep(func_name,',',''); %remove ,
        else
            func_name = 'nonlinearSys';
        end
        %generate parent object
        obj@contDynamics(func_name,ones(varargin{1},1),ones(varargin{2},1),1); %instantiate parent class
        %5 inputs
        if nargin==4
            obj.mFile = varargin{3};
            options = varargin{4};
            % link jacobian, hessian, and third order tensor files
            str = ['obj.jacobian = @jacobian_',func_name,';'];
            eval(str);
            str = ['obj.hessian = @hessianTensor_',func_name,';'];
            eval(str);
            str = ['obj.thirdOrderTensor = @thirdOrderTensor_',func_name,';'];
            eval(str);
            for i = 4:options.tensorOrder
                str = sprintf('obj.tensors{%i} = @tensor%i_%s;',i-3,i,func_name);
                eval(str);
            end
            %symbolic computation of the derivatives of the nonlinear system 
            derivatives(obj,options);
        end
    end
    
end

methods(Access = protected)
    % Override copyElement method:
    function cp = copyElement(obj)
        cp = nonlinearSys(length(obj.stateIDs), length(obj.inputIDs), str2func(obj.name));
        cp.mFile = obj.mFile;
        cp.jacobian = obj.jacobian;
        cp.linError = obj.linError;
    end
end
end


%------------- END OF CODE --------------