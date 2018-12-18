classdef nonlinearSysDT < contDynamics
% nonlinearSysDT class (nonlinear system with discrete time)
%
% Syntax:  
%    object constructor: obj = nonlinearSysDT(dim,nrOfInputs,dynFile,options)
%
% Inputs:
%    dim - system dimension
%    nrOfInputs - number of inputs
%    dynFile - handle of the dynamics file
%    options - struct containing the algorithm options
%
% Outputs:
%    obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    mFile = [];
    jacobian = [];
    hessian = [];
    thirdOrderTensor = [];
    linError = [];
end

methods
    %class constructor
    function obj = nonlinearSysDT(varargin)
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
            % link jacobian, hessian, and third order tensor files
            str = ['obj.jacobian = @jacobian_',func_name,';'];
            eval(str);
            str = ['obj.hessian = @hessianTensor_',func_name,';'];
            eval(str);
            str = ['obj.thirdOrderTensor = @thirdOrderTensor_',func_name,';'];
            eval(str);
            %symbolic computation of the derivatives of the nonlinear system 
            derivatives(obj,varargin{4});
        end
    end
   
         
    %methods in seperate files 
    [Rnext, options] = reach(obj, options)
    [obj, t, x, index] = simulate(obj, opt, tstart, tfinal, x0, y0, options)
    res = simulate_random(obj, options, runs, fractionVertices, fractionInputVertices)
    
    %display functions
    display(obj)
end
end

%------------- END OF CODE --------------