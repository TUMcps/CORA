classdef contDynamics < matlab.mixin.Copyable
% contDynamics class - basic class for continuous dynamics
%
% Syntax:  
%    object constructor: Obj = contDynamics(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    name - name of the continuous dynamics: char array
%    stateIDs - stateIDs: int array
%    inputIDs - inputIDs: int array
%    outputIDs - outputIDs: int array
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
% Written:      02-May-2007 
% Last update:  18-March-2016
%               04-August-2016 (changed to new OO format)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    name = [];
    stateIDs = [];
    inputIDs = [];
    outputIDs = [];
    dim = [];
    nrOfInputs = [];
end
    
methods
    %class constructor
    function obj = contDynamics(varargin)
        %1 input
        if nargin==1
            obj.name = varargin{1};
        %4 inputs
        elseif nargin==4
            obj.name = varargin{1};
            obj.stateIDs = varargin{2};
            obj.inputIDs = varargin{3};
            obj.outputIDs = varargin{4};
            obj.dim = length(obj.stateIDs);
            obj.nrOfInputs = length(obj.inputIDs);
        end
    end


%         %link jacobian and hessian files
%         str = ['obj.jacobian = @jacobian_',func2str(obj.dynFile),';'];
%         eval(str);
%         str = ['obj.hessian = @hessianTensor_',func2str(obj.dynFile),';'];
%         eval(str);
%         str = ['obj.hessianAbs = @hessianTensor_abs_',func2str(obj.dynFile),';'];
%         eval(str);
%         str = ['obj.thirdOrderTensor = @thirdOrderTensor_',func2str(obj.dynFile),';'];
%         eval(str);
%             
%             %compute derivatives
%             symbolicDerivation(obj,options);
%             try 
%                 if strcmp(options.category,'powerSystem')
%                     obj.other.outputBuses = options.outputBuses;
%                     obj.other.inputBuses = options.inputBuses;
%                     obj.other.buses = options.buses;
%                     obj.other.generatorBuses = options.generatorBuses;
%                     obj.other.faultBus = options.faultBus;
%                     obj.other.slackBusFlag = options.slackBusFlag;
%                     obj.other.Vgen = options.Vgen;
%                 end
%             catch
%             end
%         end
        
         
    %methods in seperate files 
    derivatives(varargin) % computes multivariate derivatives (jacobians, hessians, etc.) of nonlinear systems in a symbolic way
    [dim]=dimension(obj) % display - Returns the dimension of the system
    display(obj)  % display functions
    [Rcont,Rcont_tp,Rcont_y,Rout,Rout_tp] = reach(obj,options) % computes the reachable continuous set for the entire time horizon of a continuous system
    res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices, inputChanges) % generates random simulation of a dynamic system
    [X_compl,xTraj_compl] = simulate_rrt(obj, options, nrOfSamples, Rcont, extremePointSampling, stretchingFactor) % simulate_rrt - simulates a system using rapidly exploring random trees
    [x,u,y,p,dx,du,dy] = symVariables(varargin) % creates symbolic variables for symbolic computations
end

end

%------------- END OF CODE --------------