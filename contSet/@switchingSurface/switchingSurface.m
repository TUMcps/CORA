classdef switchingSurface
% switchingSurface class 
%
% Syntax:  
%    object constructor: Obj = switchingSurface(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    mFile - handle to an implicit function h(x) determining the switching
%    surface
%    timeOrder - polynomial order of time model for switching time
%    stateOrder - polynomial order to describe the state trajectory for
%    determining where the switching surface is hit
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    mFile = [];
    timeOrder = 2;
    stateOrder = 2;
    dim = [];
    timeModelFile = [];
    jacobianFile = [];
    nrOfParameters = [];
    y_sym = []; %store symbolic values for efficiency
    y = [];
    taylorTensors = [];
end
    
methods
    %class constructor
    function obj = switchingSurface(input1,input2,input3,input4,options)
        if nargin==1
            if isa(input1,'switchingSurface') %dummy constructor when halfspace function is called
                obj = input1;     
            end
        %one input
        elseif nargin==5
            %set handle to implicit switching surface function
            obj.mFile = input1;
            %set time order
            obj.timeOrder = input2;
            %set state order
            obj.stateOrder = input3;
            %set system dimension
            obj.dim = input4;
            %create intersection time model
            obj = intersectionTimeModel(obj,options);
            
            %link intersectionTimeModel and jacobian files
            str = ['obj.timeModelFile = @intersectionTimeModel_',func2str(obj.mFile),';'];
            eval(str);
            str = ['obj.jacobianFile = @jacobian_',func2str(obj.mFile),';'];
            eval(str);
        end
    end
         
    %methods in seperate files 
    obj = intersectionTimeModel(obj,options)
    obj = intersectionTimeParameters(obj,A,u,x_c,x_0_set)
    %obj = computeTaylorTensors(obj,A,u,x_c)
    [obj,L,Q,C,Ltilde,Qtilde] = computeTaylorTensors(obj,A,u,x_c)
    X_0 = plotPointCloud(obj,A,u,x_c,R_0,dims)
    
        
    %display functions
    display(obj)

end
end

%------------- END OF CODE --------------