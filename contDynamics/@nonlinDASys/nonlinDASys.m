classdef nonlinDASys < contDynamics
% nonlinDASys class (nonlinear differential algebraic system)
%
% Syntax:  
%    object constructor: Obj = nonlinDASys(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    dim - system dimension
%    nrOfConstraints - number of constraints
%    nrOfInputs - number of inputs
%    dynFile - handle of the dynamics file
%    conFile - handle of the constraint file
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
% Written:      27-October-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    nrOfConstraints = 0;
    dynFile = [];
    conFile = [];
    jacobian = [];
    hessian = [];
    hessianAbs = [];
    thirdOrderTensor = [];
    linError = [];
    other = [];
end
    
methods
    %class constructor
    function obj = nonlinDASys(dim,nrOfConstraints,nrOfInputs,dynFile,conFile,options)
        obj@contDynamics(func2str(dynFile),ones(dim,1),ones(nrOfInputs,1),1); %instantiate parent class
        %six inputs
        if nargin==6
            obj.nrOfConstraints=nrOfConstraints;
            obj.dynFile=dynFile;
            obj.conFile=conFile;
            
            %link jacobian and hessian files
            str = ['obj.jacobian = @jacobian_',func2str(obj.dynFile),';'];
            eval(str); 
            str = ['obj.hessian = @hessianTensor_',func2str(obj.dynFile),';'];
            eval(str);
            str = ['obj.hessianAbs = @hessianTensor_abs_',func2str(obj.dynFile),';'];
            eval(str);
            str = ['obj.thirdOrderTensor = @thirdOrderTensor_',func2str(obj.dynFile),';'];
            eval(str);
            
            %compute derivatives
            %symbolicDerivation(obj,options);
            derivatives(obj,options);
            try 
                if strcmp(options.category,'powerSystem')
                    obj.other.outputBuses = options.outputBuses;
                    obj.other.inputBuses = options.inputBuses;
                    obj.other.buses = options.buses;
                    obj.other.generatorBuses = options.generatorBuses;
                    obj.other.faultBus = options.faultBus;
                    obj.other.slackBusFlag = options.slackBusFlag;
                    obj.other.Vgen = options.Vgen;
                end
            catch
            end
        end
        
    end
         
    %methods in seperate files 
    [obj, Rfirst, Rfirst_y, options] = initReach(obj, Rinit, Rinit_y, options)
    [Rnext, options] = post(obj, R, options)
    [obj, t, x, index] = simulate(obj, opt, tstart, tfinal, x0, y0, options)
    [X_full,xTraj] = simulate_rrt(obj,options)
    [X_full,xTraj] = simulate_rrt_det(obj,options)
    handle = getfcn(obj, options)
    
    %display functions
    display(obj)
end
end

%------------- END OF CODE --------------