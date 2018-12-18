classdef nonlinParamSys < contDynamics
% nonlinParamSys class (nonlinear parametric system; parameters can be 
% constant or vary over time)
%
% Syntax:  
%    object constructor: Obj = nonlinParamSys(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    A - system matrix
%    B - input matrix
%    stepSize - time increment
%    taylorTerms - number of considered Taylor terms
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
% Written:      23-September-2010
% Last update:  27-October-2011
%               16-August-2016
%               02-June-2017
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
    constParam = 1; %flag if parameters are constant or time-varying 
    derivative = [];
    linError = [];
end
    
methods
    %class constructor
    function obj = nonlinParamSys(varargin)%(dim,nrOfInputs,nrOfParam,mFile,allowedError,options)
        % obtain function name
        if nargin>=3
            func_name = func2str(varargin{4});
            func_name = strrep(func_name,'@',''); %remove @
            func_name = strrep(func_name,'(',''); %remove (
            func_name = strrep(func_name,')',''); %remove )
            func_name = strrep(func_name,',',''); %remove ,
        else
            func_name = 'nonlinearParamSys';
        end
        %generate parent object
        obj@contDynamics(func_name,ones(varargin{1},1),ones(varargin{2},1),1); %instantiate parent class
        %5 inputs
        if nargin==5
            obj.nrOfParam = varargin{3};
            obj.mFile = varargin{4};
            % link jacobian, hessian, and third order tensor files
            str = ['obj.jacobian = @jacobian_',func_name,';'];
            eval(str);
            str = ['obj.jacobian_freeParam = @jacobian_freeParam_',func_name,';'];
            eval(str);
            str = ['obj.hessian = @hessianTensor_',func_name,';'];
            eval(str);
            str = ['obj.thirdOrderTensor = @thirdOrderTensor_',func_name,';'];
            eval(str);
            str = ['obj.parametricDynamicFile = @parametricDynamicFile_',func_name,';'];
            eval(str);
            %symbolic computation of the jacobians of the nonlinear system 
            %obj = symbolicDerivation(obj,varargin{5});
            derivatives(obj,varargin{5});
        end
    end
         
    %methods in seperate files 
    [obj,Rfirst,options] = initReach(obj, Rinit, options)
    [obj,t,x,index] = simulate(obj,opt,tstart,tfinal,x0,options)
    handle = getfcn(obj,options)
    
    %display functions
    display(obj)
end
end

%------------- END OF CODE --------------