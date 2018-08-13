classdef linParamSys < contDynamics
% linParamSys class (linear parametric system; parameters are constant over time)
%
% Syntax:  
%    object constructor: Obj = linParamSys(varargin)
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
% Last update:  01-November-2017 (add possibility to change between constant and varying parameters)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    A = 1;
    B = 0;
    %stepSize = 0.01/max(abs(eig(A.center)));
    constParam = 1; %flag if parameters are constant or time-varying
    stepSize = 1;
    taylorTerms = [];
    mappingMatrixSet = [];
    power = [];
    E = [];
    F = [];
    inputF = [];
    inputCorr = [];
    Rinput = [];
    Rtrans = [];
    RV = [];
    sampleMatrix = [];
end
    
methods
    %class constructor
    function obj = linParamSys(A,B,stepSize,taylorTerms,paramType)
        obj@contDynamics('linParamSysDefault',ones(A.dim,1),1,1); %instantiate parent class
        %one input
        if nargin==1
            obj.A = A;
        %two inputs
        elseif nargin==2
            obj.A = A;
            obj.B = B;
        %three inputs
        elseif nargin==3
            obj.A = A;
            obj.B = B;
            obj.stepSize = stepSize;
        %four inputs
        elseif nargin==4
            obj.A = A;
            obj.B = B;
            obj.stepSize = stepSize;
            obj.taylorTerms = taylorTerms;
        %five inputs
        elseif nargin==5
            obj.A = A;
            obj.B = B;
            obj.stepSize = stepSize;
            obj.taylorTerms = taylorTerms;
            if strcmp(paramType,'varParam')
                obj.constParam = 0;
            elseif strcmp(paramType,'constParam')
                obj.constParam = 1;
            end
        end
        
        tmp = randomSampling(obj.A,1);
        obj.sampleMatrix.A = tmp{1};
    end
         
    %methods in seperate files 
    [obj,Rfirst,options] = initReach(obj, Rinit, options)
    [Rnext,options] = post(obj,R,options)
    [obj] = preReach(obj,options)
    [Rfirst] = coreReach(obj,Rinit,options)
    [Rnext,IH] = postReach(obj,Rinit,R_tp,c)
    [obj,t,x,index] = simulate(obj,opt,tstart,tfinal,x0,options)
    handle = getfcn(obj,options)
    
    %display functions
    plot(varargin)
    display(obj)
end
end

%------------- END OF CODE --------------